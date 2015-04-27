
#ifdef HAVE_PETSCEXT

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsnes.h>
#include <petscext.h>
#include <petscext_pc.h>

#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include "petsc-private/kspimpl.h"   /*I "petscksp.h" I*/
#else
  #include "private/kspimpl.h"   /*I "petscksp.h" I*/
#endif

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "Solvers/SLE/SLE.h" /* to give the AugLagStokes_SLE type */
#include "Solvers/KSPSolvers/KSPSolvers.h"

#include "BSSCR.h"
#include "createK2.h"
#include "writeMatVec.h"
#include "ksp_pressure_nullspace.h"

#define KSPBSSCR         "bsscr"

EXTERN_C_BEGIN
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT KSPCreate_BSSCR(KSP);
EXTERN_C_END

PetscErrorCode BSSCR_DRIVER_auglag( KSP ksp, Mat stokes_A, Vec stokes_x, Vec stokes_b, Mat approxS, KSP ksp_K,
                                    MatStokesBlockScaling BA, PetscTruth sym, KSP_BSSCR * bsscrp_self );
PetscErrorCode BSSCR_DRIVER_flex( KSP ksp, Mat stokes_A, Vec stokes_x, Vec stokes_b, Mat approxS, KSP ksp_K,
                                  MatStokesBlockScaling BA, PetscTruth sym, KSP_BSSCR * bsscrp_self );

/*************************************************************************************************/
/****** START Helper Functions *******************************************************************/
#define BSSCR_GetPetscMatrix( matrix ) ( (Mat)(matrix) )
/****** END Helper Functions *********************************************************************/
/*************************************************************************************************/
typedef struct {
    void *ctx;
    KSP_BSSCR * bsscr;
} BSSCR_KSPConverged_Ctx;


/*************************************************************************************************/
/*************************************************************************************************/

#undef __FUNCT__  
#define __FUNCT__ "BSSCR_KSPSetConvergenceMinIts"
PetscErrorCode BSSCR_KSPSetConvergenceMinIts(KSP ksp, PetscInt n, KSP_BSSCR * bsscr)
{
      BSSCR_KSPConverged_Ctx *ctx;
      PetscErrorCode ierr;

      PetscFunctionBegin;

      bsscr->min_it = n; /* set minimum its */

#if(PETSC_VERSION_MAJOR == 3)
      ierr = Stg_PetscNew(BSSCR_KSPConverged_Ctx,&ctx);CHKERRQ(ierr);
      ierr = KSPDefaultConvergedCreate(&ctx->ctx);CHKERRQ(ierr);
      ctx->bsscr=bsscr;
      ierr = KSPSetConvergenceTest(ksp,BSSCR_KSPConverged,ctx,BSSCR_KSPConverged_Destroy);CHKERRQ(ierr);
#else
      ierr = KSPSetConvergenceTest(ksp,BSSCR_KSPConverged,(void*)bsscr);CHKERRQ(ierr);
#endif

      PetscFunctionReturn(0);   
}

#undef __FUNCT__  
#define __FUNCT__ "BSSCR_KSPConverged"
PetscErrorCode BSSCR_KSPConverged(KSP ksp,PetscInt n,PetscReal rnorm,KSPConvergedReason *reason,void *cctx)
{
  PetscErrorCode           ierr;
#if(PETSC_VERSION_MAJOR == 3)
  BSSCR_KSPConverged_Ctx  *ctx = (BSSCR_KSPConverged_Ctx *)cctx;
  KSP_BSSCR               *bsscr = ctx->bsscr;
#else
  KSP_BSSCR               *bsscr = (KSP_BSSCR*)cctx;
#endif


  PetscFunctionBegin;
#if(PETSC_VERSION_MAJOR == 3)
  ierr = KSPDefaultConverged(ksp,n,rnorm,reason,ctx->ctx);CHKERRQ(ierr);
#else
  ierr = KSPDefaultConverged(ksp,n,rnorm,reason,cctx);CHKERRQ(ierr);
#endif
  if (*reason) {
    ierr = PetscInfo2(ksp,"default convergence test KSP iterations=%D, rnorm=%G\n",n,rnorm);CHKERRQ(ierr);
  }
  if(ksp->its < bsscr->min_it){
      ksp->reason          = KSP_CONVERGED_ITERATING;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BSSCR_KSPConverged_Destroy"
PetscErrorCode BSSCR_KSPConverged_Destroy(void *cctx)
{
  BSSCR_KSPConverged_Ctx *ctx = (BSSCR_KSPConverged_Ctx *)cctx;
  PetscErrorCode           ierr;

  PetscFunctionBegin;
  ierr = KSPDefaultConvergedDestroy(ctx->ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSPRegisterBSSCR"
PetscErrorCode PETSCKSP_DLLEXPORT KSPRegisterBSSCR(const char path[])
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = Stg_KSPRegister(KSPBSSCR, path, "KSPCreate_BSSCR", KSPCreate_BSSCR );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "KSPSolve_BSSCR"
PetscErrorCode  KSPSolve_BSSCR(KSP ksp)
{
    Mat            Amat,Pmat; /* Stokes Matrix and it's Preconditioner matrix: Both 2x2 PetscExt block matrices */
    Vec            B, X; /* rhs and solution vectors */
    MatStructure   pflag;
    PetscErrorCode ierr;
    KSP_BSSCR *    bsscr;    
    Stokes_SLE *    SLE;
    PETScMGSolver * MG;
    Mat K,D,ApproxS;
    KSP A11_ksp; /* Velocity Solve KSP */
    MatStokesBlockScaling BA;
    PetscTruth flg, sym, augment;
    
    PetscFunctionBegin;
    PetscPrintf( PETSC_COMM_WORLD, "**** BSSCR -- Block Stokes Schur Compliment Reduction Solver **** \n");
    /** Get the stokes Block matrix and its preconditioner matrix */
    ierr = Stg_PCGetOperators(ksp->pc,&Amat,&Pmat,&pflag);CHKERRQ(ierr);
    /** In Petsc proper, KSP's ksp->data is usually set in KSPCreate_XXX function. 
        Here it is set in the _StokesBlockKSPInterface_Solve function instead so that we can ensure that the solver
        has everything it needs */
    bsscr         = (KSP_BSSCR*)ksp->data;
    MG            = (PETScMGSolver*)bsscr->mg;
    SLE           = (Stokes_SLE*)bsscr->st_sle;
    X             = ksp->vec_sol;
    B             = ksp->vec_rhs;

    if( bsscr->do_scaling ){
        (*bsscr->scale)(ksp); /* scales everything including the UW preconditioner */
        BA =  bsscr->BA;
    }

    if( bsscr->k2type != PETSC_NULL && bsscr->buildK2 ){
        (*bsscr->buildK2)(ksp); /* building K2 from scaled version of stokes operators: K2 lives on bsscr struct = ksp->data */  
    }

    /* get sub matrix / vector objects */
    MatNestGetSubMat( Amat, 0,0, &K );
    /* Underworld preconditioner matrix*/
    ApproxS = PETSC_NULL;
    if( ((StokesBlockKSPInterface*)SLE->solver)->preconditioner ) { /* SLE->solver->st_sle == SLE here, by the way */
        StiffnessMatrix *preconditioner;
        preconditioner = ((StokesBlockKSPInterface*)SLE->solver)->preconditioner;
        ApproxS = BSSCR_GetPetscMatrix( preconditioner->matrix ); 
    }

    sym = bsscr->DIsSym;

    KSPCreate(PETSC_COMM_WORLD, &A11_ksp);
    Stg_KSPSetOperators(A11_ksp, K, K, DIFFERENT_NONZERO_PATTERN);

    MatNestGetSubMat( Amat, 1,0, &D );if(!D){ PetscPrintf( PETSC_COMM_WORLD, "D does not exist but should!!\n"); exit(1); }

    /**********************************************************/
    /******* SOLVE!! ******************************************/
    /**********************************************************/
    flg = PETSC_FALSE;
    augment = PETSC_TRUE;
    PetscOptionsGetTruth(PETSC_NULL, "-augmented_lagrangian", &augment, &flg);
    if(augment){
        BSSCR_DRIVER_auglag( ksp, Amat, X, B, ApproxS, A11_ksp, BA, sym, bsscr ); }
    else{
        BSSCR_DRIVER_flex( ksp, Amat, X, B, ApproxS, A11_ksp, BA, sym, bsscr ); }

    /**********************************************************/
    /***** END SOLVE!! ****************************************/
    /**********************************************************/
    if( bsscr->do_scaling ){ 
        (*bsscr->unscale)(ksp);  }
    if( bsscr->k2type != PETSC_NULL && bsscr->K2 != PETSC_NULL ){
        if(bsscr->k2type != K2_SLE){/* don't destroy here, as in this case, K2 is just pointing to an existing matrix on the SLE */
            Stg_MatDestroy(&bsscr->K2 );
        }
        bsscr->K2built = PETSC_FALSE;
        bsscr->K2 = PETSC_NULL;  }

    ksp->reason = KSP_CONVERGED_RTOL;

    PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPDestroy_BSSCR" 
PetscErrorCode KSPDestroy_BSSCR(KSP ksp)
{
    KSP_BSSCR     *bsscr = (KSP_BSSCR *)ksp->data;
    Mat            K2=bsscr->K2;
    Vec            t,v;
    MatStokesBlockScaling BA=bsscr->BA;
    PetscErrorCode ierr;

    PetscFunctionBegin;

    t=bsscr->t;
    v=bsscr->v;
    if ( v ){ Stg_VecDestroy(&v); }
    if ( t ){ Stg_VecDestroy(&t); }
    if( K2 ){ Stg_MatDestroy(&K2 ); }/* shouldn't need this now */
    if(BA) BSSCR_MatStokesBlockScalingDestroy( BA );
    ierr = PetscFree(ksp->data);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
static const char *K2Types[] = {"NULL","DGMGD","GMG","GG","SLE","K2Type","K2_",0};
#undef __FUNCT__  
#define __FUNCT__ "KSPSetFromOptions_BSSCR"
PetscErrorCode KSPSetFromOptions_BSSCR(KSP ksp)
{

    PetscTruth  flg;
    KSP_BSSCR  *bsscr = (KSP_BSSCR *)ksp->data;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = PetscOptionsHead("KSP BSSCR options");CHKERRQ(ierr);
    /* if this ksp has a prefix "XXX_" it will be automatically added to the options. e.g. -ksp_test -> -XXX_ksp_test */
    /* ierr = PetscOptionsTruth("-ksp_test","Test KSP flag","nil",PETSC_FALSE,&test,PETSC_NULL);CHKERRQ(ierr); */
    /* if(test){ PetscPrintf( PETSC_COMM_WORLD,  "\n\n-----  test flag set  ------\n\n"); } */
    ierr = PetscOptionsEnum("-ksp_k2_type","Augmented Lagrangian matrix type","",K2Types, bsscr->k2type,(PetscEnum*)&bsscr->k2type,&flg);CHKERRQ(ierr);
    //if(flg){  PetscPrintf( PETSC_COMM_WORLD,  "-----  k2 type is  ------\n"); }
    ierr = PetscOptionsTail();CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPView_BSSCR" 
PetscErrorCode KSPView_BSSCR(KSP ksp,PetscViewer viewer)
{
    PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "KSPSetUp_BSSCR" 
PetscErrorCode KSPSetUp_BSSCR(KSP ksp)
{
    Vec            B,X;
    KSP_BSSCR  *bsscr = (KSP_BSSCR *)ksp->data;
    Mat K,G, ApproxS;
    Stokes_SLE*  stokesSLE  = (Stokes_SLE*)bsscr->st_sle;
    PetscTruth sym,flg,ismumps,augment,scale,konly,found,conp,checkerp;
    char name[PETSC_MAX_PATH_LEN];

    PetscFunctionBegin;

    BSSCR_PetscExtStokesSolversInitialize();

    ApproxS = bsscr->preconditioner->matrix;

    K = stokesSLE->kStiffMat->matrix;
    G = stokesSLE->gStiffMat->matrix;

    bsscr->K2=PETSC_NULL;
    X             = ksp->vec_sol;
    B             = ksp->vec_rhs;

    BSSCR_MatStokesBlockScalingCreate( &(bsscr->BA) );/* allocate memory for scaling struct */
    found = PETSC_FALSE;
    augment = PETSC_TRUE;
    PetscOptionsGetTruth(PETSC_NULL, "-augmented_lagrangian", &augment, &found);
    //PetscOptionsGetString( PETSC_NULL, "-no_augmented_lagrangian", name, PETSC_MAX_PATH_LEN-1, &flg );
    if(augment){
        bsscr->buildK2 = bsscr_buildK2; }
    else {
        bsscr->buildK2 = PETSC_NULL; }

    sym = bsscr->DIsSym;

    /***************************************************************************************************************/
    /** Do scaling *************************************************************************************************/
    /***************************************************************************************************************/
    found = PETSC_FALSE;
    scale = PETSC_FALSE;/* scaling off by default */
    PetscOptionsGetTruth(PETSC_NULL, "-rescale_equations", &scale, &found);
    /* if -rescale_equations option given then found and scale are PETSC_TRUE */
    /* keeping the -no_scale option around for time being if above option not given */
    //if(!found){
    //    flg = PETSC_FALSE;
    //    scale = PETSC_TRUE;/* scaling on by default */
    //    PetscOptionsGetString( PETSC_NULL, "-no_scale", name, PETSC_MAX_PATH_LEN-1, &flg );
    //    if(flg){ scale = PETSC_FALSE; }
    //}
    //Stg_PetscObjectTypeCompare((PetscObject)K, MATAIJMUMPS, &ismumps);
    Stg_PetscObjectTypeCompare((PetscObject)K, "mataijmumps", &ismumps);/** older versions of petsc have this */
    if(ismumps && scale){
        PetscPrintf( PETSC_COMM_WORLD, "\t* Not applying scaling to matrices as MatGetRowMax operation not defined for MATAIJMUMPS matrix type \n");
        scale = PETSC_FALSE;
    }
    if( scale ) {
        bsscr->scale       = KSPScale_BSSCR;
        bsscr->unscale     = KSPUnscale_BSSCR;
        bsscr->do_scaling  = PETSC_TRUE; 
        bsscr->scaletype   = KONLY;
        konly = PETSC_TRUE;
        found = PETSC_FALSE;
        PetscOptionsGetTruth(PETSC_NULL, "-k_scale_only", &konly, &found);
        //PetscOptionsGetString( PETSC_NULL, "-k_scale", name, PETSC_MAX_PATH_LEN-1, &flg2 );/* scaling on K only*/
        if(!konly){ bsscr->scaletype = DEFAULT; }
    }
    //else {
    //    PetscPrintf( PETSC_COMM_WORLD, "\t* Not applying any scaling to matrices \n");
    //}
    /***************************************************************************************************************/
    /**  Set up functions for building Pressure Null Space Vectors *************************************************/
    /***************************************************************************************************************/
    found = PETSC_FALSE;
    checkerp = PETSC_FALSE;
    //PetscOptionsGetString( PETSC_NULL, "-build_cb_pressure_nullspace", name, PETSC_MAX_PATH_LEN-1, &flg );
    PetscOptionsGetTruth(PETSC_NULL, "-remove_checkerboard_pressure_null_space", &checkerp, &found );
    if(checkerp){
        bsscr->check_cb_pressureNS    = PETSC_TRUE;
        bsscr->check_pressureNS       = PETSC_TRUE;
        bsscr->buildPNS               = KSPBuildPressure_CB_Nullspace_BSSCR;
    }
    //flg = PETSC_FALSE;
    //PetscOptionsGetString( PETSC_NULL, "-build_const_pressure_nullspace", name, PETSC_MAX_PATH_LEN-1, &flg );
    found = PETSC_FALSE;
    conp = PETSC_FALSE;
    PetscOptionsGetTruth(PETSC_NULL, "-remove_constant_pressure_null_space", &conp, &found );
    if(conp){
        bsscr->check_const_pressureNS = PETSC_TRUE;
        bsscr->check_pressureNS       = PETSC_TRUE;
        bsscr->buildPNS               = KSPBuildPressure_Const_Nullspace_BSSCR;
    }
    /***************************************************************************************************************/
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "KSPCreate_BSSCR"
PetscErrorCode PETSCKSP_DLLEXPORT KSPCreate_BSSCR(KSP ksp)
{
    KSP_BSSCR  *bsscr;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = Stg_PetscNew(KSP_BSSCR,&bsscr);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory(ksp,sizeof(KSP_BSSCR));CHKERRQ(ierr);
    //ierr = PetscNewLog(ksp,KSP_BSSCR,&bsscr);CHKERRQ(ierr);
    ksp->data                              = (void*)bsscr;
    
    #if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2 ) )
    ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,0);CHKERRQ(ierr);
    ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_LEFT,1);CHKERRQ(ierr);
    ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NATURAL,PC_LEFT,0);CHKERRQ(ierr);
    ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_LEFT,1);CHKERRQ(ierr);
    #endif
    /*
       Sets the functions that are associated with this data structure 
       (in C++ this is the same as defining virtual functions)
    */
    ksp->ops->setup                = KSPSetUp_BSSCR;
    ksp->ops->solve                = KSPSolve_BSSCR;
    ksp->ops->destroy              = KSPDestroy_BSSCR;
    ksp->ops->view                 = KSPView_BSSCR;
    ksp->ops->setfromoptions       = KSPSetFromOptions_BSSCR;
    ksp->ops->buildsolution        = KSPDefaultBuildSolution;
    ksp->ops->buildresidual        = KSPDefaultBuildResidual;

//    bsscr->k2type=K2_GMG;
    bsscr->k2type      = 0;
    bsscr->do_scaling  = PETSC_FALSE;
    bsscr->scaled      = PETSC_FALSE;
    bsscr->K2built     = PETSC_FALSE;
    bsscr->check_cb_pressureNS     = PETSC_FALSE;/* checker board nullspace */
    bsscr->check_const_pressureNS  = PETSC_FALSE;/* constant nullspace */
    bsscr->check_pressureNS        = PETSC_FALSE;/* is true if either of above two are true */
    bsscr->t           = NULL;/* null space vectors for pressure */
    bsscr->v           = NULL;
    bsscr->nstol       = 1e-7;/* null space detection tolerance */
    PetscFunctionReturn(0);
}
EXTERN_C_END

#endif
