#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscis.h> 

#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include <petsc-private/kspimpl.h>
#else
  #include <private/kspimpl.h>
#endif


//#include "common-driver-utils.h"

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>
#include "Solvers/SLE/SLE.h" /* to give the AugLagStokes_SLE type */
#include "Solvers/KSPSolvers/KSPSolvers.h"
#include "BSSCR.h"
#include "stokes_block_scaling.h"
#include "stokes_mvblock_scaling.h"
#include "writeMatVec.h"
#include "createK2.h"
#include "mg.h"
#include "summary.h"
#include "ksp_pressure_nullspace.h"

#define BSSCR_GetPetscMatrix( matrix ) ( (Mat)(matrix) )
#define BSSCR_GetPetscVector( vector ) ( (Vec)(vector) )

typedef struct {
    PetscTruth initialrtol;    /* default relative residual decrease is computing from initial residual, not rhs */
    PetscTruth mininitialrtol; /* default relative residual decrease is computing from min of initial residual and rhs */
    Vec        work;
    PetscReal  pointwise_max;
} KSPPWConvergedCtx;

PetscErrorCode BSSCR_KSPNormInfMonitor( KSP ksp, PetscInt iteration, PetscReal residualNorm, void *dummy);
PetscErrorCode BSSCR_KSPNormInfToNorm2Monitor( KSP ksp, PetscInt iteration, PetscReal residualNorm, void *dummy);

PetscErrorCode BSSCR_KSPPWConvergedCreate(void **ctx);
PetscErrorCode BSSCR_KSPPWConvergedDestroy(void *ctx);
PetscErrorCode BSSCR_KSPNormInfConverged(KSP ksp,PetscInt n,PetscReal rnorm,KSPConvergedReason *reason,void *ctx);
PetscErrorCode BSSCR_KSPNormInfToNorm2Monitor(KSP ksp,PetscInt n,PetscReal rnorm, void *dummy);
PetscErrorCode BSSCR_KSPNormInfMonitor(KSP ksp,PetscInt n,PetscReal rnorm, void *dummy);
PetscErrorCode BSSCR_KSPNorm2RawMonitor(KSP ksp,PetscInt n,PetscReal rnorm, void *dummy);
PetscErrorCode BSSCR_KSPSetNormInfConvergenceTest(KSP ksp);


#undef __FUNCT__  
#define __FUNCT__ "BSSCR_DRIVER_auglag"
PetscErrorCode BSSCR_DRIVER_auglag( KSP ksp, Mat stokes_A, Vec stokes_x, Vec stokes_b, Mat approxS,
                                          MatStokesBlockScaling BA, PetscTruth sym, KSP_BSSCR * bsscrp_self )
{
    AugLagStokes_SLE *    stokesSLE = (AugLagStokes_SLE*)bsscrp_self->st_sle;
    PetscTruth uzawastyle, KisJustK=PETSC_TRUE, restorek, change_A11rhspresolve;
    PetscTruth usePreviousGuess, useNormInfStoppingConditions, useNormInfMonitor, found, forcecorrection;
    PetscTruth change_backsolve;
    PetscErrorCode ierr;
    PetscInt monitor_index,max_it,min_it;
    KSP ksp_inner, ksp_S, ksp_new_inner;
    PC pc_S, pcInner;
    Mat K,G,D,C, S, K2;// Korig;
    Vec u,p,f,f2,f3=0,h, h_hat,t;
    IS  isr[2];
    MGContext mgCtx;
    double mgSetupTime, scrSolveTime, a11SingleSolveTime, penaltyNumber, hFactor;
    static int been_here = 0;  /* Ha Ha Ha !! */

    char name[PETSC_MAX_PATH_LEN];
    char matname[PETSC_MAX_PATH_LEN];
    char suffix[PETSC_MAX_PATH_LEN];
    char str[PETSC_MAX_PATH_LEN];
    PetscTruth flg, extractMats;
    
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    if( bsscrp_self->st_sle->context->loadFromCheckPoint ){
        been_here=1;
    }
    /* get sub matrix / vector objects */
    /* note that here, the matrix D should always exist. It is set up in  _StokesBlockKSPInterface_Solve in StokesBlockKSPInterface.c */
    /* now extract K,G etc from a MatNest object */
    /* {VecGetSubVector(y,bA->isglobal.row[i],&by[i]);} */

    MatNestGetSubMat( stokes_A, 0,0, &K );
    MatNestGetSubMat( stokes_A, 0,1, &G );
    MatNestGetSubMat( stokes_A, 1,0, &D );if(!D){ PetscPrintf( PETSC_COMM_WORLD, "D does not exist but should!!\n"); exit(1); }
    MatNestGetSubMat( stokes_A, 1,1, &C );  
    VecNestGetSubVec( stokes_x, 0, &u );
    VecNestGetSubVec( stokes_x, 1, &p );
    VecNestGetSubVec( stokes_b, 0, &f );
    VecNestGetSubVec( stokes_b, 1, &h );
    
    PetscPrintf( PETSC_COMM_WORLD,  "\n\n----------  AUGMENTED LAGRANGIAN K2 METHOD ---------\n\n" );
    PetscPrintf( PETSC_COMM_WORLD,      "----------- Penalty = %f\n\n", stokesSLE->penaltyNumber );
    sprintf(suffix,"%s","x");

    PetscOptionsGetString( PETSC_NULL, "-matsuffix", suffix, PETSC_MAX_PATH_LEN-1, &extractMats );

    flg=0;
    PetscOptionsGetString( PETSC_NULL, "-matdumpdir", name, PETSC_MAX_PATH_LEN-1, &flg );
    if(flg){
        sprintf(str,"%s/",name); sprintf(matname,"K%s",suffix);
        bsscr_dirwriteMat( K, matname,str, "Writing K matrix in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"G%s",suffix);
        bsscr_dirwriteMat( G, matname,str, "Writing G matrix in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"f%s",suffix);
        bsscr_dirwriteVec( f, matname,str, "Writing f vector in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"h%s",suffix);
        bsscr_dirwriteVec( h, matname,str, "Writing h vector in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"Shat%s",suffix);
        bsscr_dirwriteMat( approxS, matname,str, "Writing Shat matrix in al Solver");
    }

    penaltyNumber = stokesSLE->penaltyNumber;
    hFactor = stokesSLE->hFactor;
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /******  GET K2   ****************************************************************************************/
    if(penaltyNumber > 1e-10 && bsscrp_self->k2type){

        flg=0;
        PetscOptionsGetString( PETSC_NULL, "-matdumpdir", name, PETSC_MAX_PATH_LEN-1, &flg );
        if(flg){
            sprintf(str,"%s/",name);   sprintf(matname,"K2%s",suffix);
            bsscr_dirwriteMat( bsscrp_self->K2, matname,str, "Writing K2 matrix in al Solver");
        }

        K2=bsscrp_self->K2;
        scrSolveTime = MPI_Wtime();
        ierr=MatAXPY(K,penaltyNumber,K2,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);/* Computes K = penaltyNumber*K2 + K */
        scrSolveTime =  MPI_Wtime() - scrSolveTime;
        PetscPrintf( PETSC_COMM_WORLD, "\n\t* K+p*K2 in time: %lf seconds\n\n", scrSolveTime);

        //ierr=MatAYPX(K2,penaltyNumber,K,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);/* Computes K2 = penaltyNumber*K2 + K */
        //Korig=K;
        //K=K2; /* we are doing K2=a*K2+K because K2s non-zero pattern is larger than Ks. The other way round is very slow due to extra memory allocation in K. */

        KisJustK=PETSC_FALSE;
        forcecorrection=PETSC_TRUE;
        PetscOptionsGetTruth( PETSC_NULL ,"-force_correction", &forcecorrection, &found );
        if(forcecorrection){
            if(bsscrp_self->f2 && forcecorrection){ 
                f2=bsscrp_self->f2;
                ierr=VecAXPY(f,penaltyNumber,f2);/*  f <- f +a*f2 */
            }else{
                switch (bsscrp_self->k2type) {
                case (K2_GG):
                {
                    VecDuplicate( u, &f3 );
                    MatMult( G, h, f3);   ierr=VecAXPY(f,penaltyNumber,f3);/*  f <- f +a*f2 */
                }
                break;
                case (K2_GMG):
                {
                    Mat M;
                    Vec Mdiag;
                    VecDuplicate( u, &f3 );
                    M = stokesSLE->mStiffMat->matrix;
                    MatGetVecs( M, &Mdiag, PETSC_NULL );
                    MatGetDiagonal( M, Mdiag );
                    VecReciprocal(Mdiag);
                    VecPointwiseMult(Mdiag,Mdiag,h);
                    MatMult( G, Mdiag, f3);   ierr=VecAXPY(f,penaltyNumber,f3);/*  f <- f +a*f2 */
                    Stg_VecDestroy(&Mdiag);
                }
                break;
                case (K2_DGMGD):
                {
                    Mat M;
                    Vec Mdiag;
                    VecDuplicate( u, &f3 );
                    M = stokesSLE->mStiffMat->matrix;
                    MatGetVecs( M, &Mdiag, PETSC_NULL );
                    MatGetDiagonal( M, Mdiag );
                    VecReciprocal(Mdiag);
                    VecPointwiseMult(Mdiag,Mdiag,h);
                    MatMult( G, Mdiag, f3);   ierr=VecAXPY(f,penaltyNumber,f3);/*  f <- f +a*f2 */
                    Stg_VecDestroy(&Mdiag);
                }
                break;
                case (K2_NULL):
                {
                    ;
                }
                break;
                case (K2_SLE):
                {
                    ;
                }
                break;
                }
            }
        }
    }

    /* Create Schur complement matrix */
    MatCreateSchurComplement(K,K,G,D,C, &S);
    //MatCreateSchurFromBlock( stokes_A, 0.0, "MatSchur_A11", &S );
    //MatAssemblyBegin( S, MAT_FINAL_ASSEMBLY );
    //MatAssemblyEnd( S, MAT_FINAL_ASSEMBLY );

    /* configure inner solver */
    //if (!ksp_K) {  PetscPrintf( PETSC_COMM_WORLD,"ksp_K cannot be NULL\n"); abort();}

    //MatSchurSetKSP( S, ksp_K );
    //MatSchurComplementSetKSP( S, ksp_K);
    //MatSchurGetKSP( S, &ksp_inner );
    MatSchurComplementGetKSP( S, &ksp_inner);
    KSPGetPC( ksp_inner, &pcInner );
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /*********    SET PREFIX FOR INNER/VELOCITY KSP    *************************************************************/
    KSPSetOptionsPrefix( ksp_inner, "A11_" );
    KSPSetFromOptions( ksp_inner );
    Stg_KSPSetOperators(ksp_inner, K, K, DIFFERENT_NONZERO_PATTERN);

    useNormInfStoppingConditions = PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-A11_use_norm_inf_stopping_condition", &useNormInfStoppingConditions, &found );
    if(useNormInfStoppingConditions) 
        BSSCR_KSPSetNormInfConvergenceTest( ksp_inner ); 

    useNormInfMonitor = PETSC_FALSE; 
    PetscOptionsGetTruth( PETSC_NULL, "-A11_ksp_norm_inf_monitor", &useNormInfMonitor, &found );
    if(useNormInfMonitor)  KSPMonitorSet( ksp_inner, BSSCR_KSPNormInfMonitor, PETSC_NULL, PETSC_NULL );

    useNormInfMonitor = PETSC_FALSE; 
    PetscOptionsGetTruth( PETSC_NULL, "-A11_ksp_norm_inf_to_norm_2_monitor", &useNormInfMonitor, &found );
    if(useNormInfMonitor)  KSPMonitorSet( ksp_inner, BSSCR_KSPNormInfToNorm2Monitor, PETSC_NULL, PETSC_NULL );
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /* If multigrid is enabled, set it now. */
    change_A11rhspresolve = PETSC_FALSE; 
    PetscOptionsGetTruth( PETSC_NULL, "-change_A11rhspresolve", &change_A11rhspresolve, &found );

    if(bsscrp_self->mg && !change_A11rhspresolve) { mgSetupTime=setupMG( bsscrp_self, ksp_inner, pcInner, K, &mgCtx ); }
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /* create right hand side */
    if(change_A11rhspresolve){
      //Stg_KSPDestroy(&ksp_inner );
      KSPCreate(PETSC_COMM_WORLD, &ksp_new_inner);
      Stg_KSPSetOperators(ksp_new_inner, K, K, DIFFERENT_NONZERO_PATTERN);
      KSPSetOptionsPrefix(ksp_new_inner, "rhsA11_");
      MatSchurComplementSetKSP( S, ksp_new_inner );/* this call destroys the ksp_inner that is already set on S */
      ksp_inner=ksp_new_inner;
      KSPGetPC( ksp_inner, &pcInner );
      KSPSetFromOptions(ksp_inner); /* make sure we are setting up our solver how we want it */
    }
    MatGetVecs( S, PETSC_NULL, &h_hat );
    //MatSchurApplyReductionToVecFromBlock( S, stokes_b, h_hat );/* A11 KSPSolve in here */
    Vec f_tmp;
    /* It may be the case that the current velocity solution might not be bad guess for f_tmp? */
    MatGetVecs( K, PETSC_NULL, &f_tmp );
    scrSolveTime = MPI_Wtime();
    KSPSolve(ksp_inner, f, f_tmp);
    scrSolveTime =  MPI_Wtime() - scrSolveTime;
    PetscPrintf( PETSC_COMM_WORLD, "\n\t*  KSPSolve for RHS setup Finished in time: %lf seconds\n\n", scrSolveTime);
    //bsscr_writeVec( t, "ts", "Writing t vector");
    MatMult(D, f_tmp, h_hat);
    //VecAXPY(h, -1, h_hat);/* h_hat = h - Gt*K^(-1)*f */
    VecAYPX(h_hat, -1.0, h); /* Computes y = x + alpha y.  h_hat -> h - Gt*K^(-1)*f*/
    Stg_VecDestroy(&f_tmp);

    if(bsscrp_self->mg && change_A11rhspresolve) {
      //Stg_KSPDestroy(&ksp_inner );
      KSPCreate(PETSC_COMM_WORLD, &ksp_new_inner);
      Stg_KSPSetOperators(ksp_new_inner, K, K, DIFFERENT_NONZERO_PATTERN);
      KSPSetOptionsPrefix( ksp_new_inner, "A11_" );
      MatSchurComplementSetKSP( S, ksp_new_inner );
      ksp_inner=ksp_new_inner;
      //MatSchurSetKSP( S, ksp_inner );
      KSPGetPC( ksp_inner, &pcInner );
      KSPSetFromOptions( ksp_inner );
      mgSetupTime=setupMG( bsscrp_self, ksp_inner, pcInner, K, &mgCtx ); 
    }
    /* create solver for S p = h_hat */
    KSPCreate( PETSC_COMM_WORLD, &ksp_S );
    KSPSetOptionsPrefix( ksp_S, "scr_");
    /* By default use the UW approxS Schur preconditioner -- same as the one used by the Uzawa solver */
    /* Note that if scaling is activated then the approxS matrix has been scaled already */
    /* so no need to rebuild in the case of scaling as we have been doing */
    if(!approxS){ PetscPrintf( PETSC_COMM_WORLD,  "WARNING approxS is NULL\n"); }
    Stg_KSPSetOperators( ksp_S, S, S, SAME_NONZERO_PATTERN );
    KSPSetType( ksp_S, "cg" );
    KSPGetPC( ksp_S, &pc_S );
    BSSCR_BSSCR_StokesCreatePCSchur2( K,G,D,C,approxS, pc_S, sym, bsscrp_self );

    flg=0;
    PetscOptionsGetString( PETSC_NULL, "-NN", name, PETSC_MAX_PATH_LEN-1, &flg );
    if(flg){
	Mat Smat, Pmat;                                                                                                 
	MatStructure mstruct;  
	Stg_PCGetOperators( pc_S, &Smat, &Pmat, &mstruct );
        sprintf(str,"%s/",name); sprintf(matname,"Pmat%s",suffix);
        bsscr_dirwriteMat( Pmat, matname,str, "Writing Pmat matrix in al Solver");
    }

    
    uzawastyle=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL, "-uzawa_style", &uzawastyle, &found );
    if(uzawastyle){
        /* now want to set up the ksp_S->pc to be of type ksp (gmres) by default to match Uzawa */
        KSP pc_ksp;
        KSPGetPC( ksp_S, &pc_S );
        PCSetType(pc_S,PCKSP);
        PCKSPGetKSP( pc_S, &pc_ksp);
        KSPSetType(pc_ksp, "gmres" );
        KSPSetOptionsPrefix( pc_ksp, "scrPCKSP_");
        KSPSetFromOptions( pc_ksp );
    }
    KSPSetFromOptions( ksp_S );
    /* Set specific monitor test */
    KSPGetTolerances( ksp_S, PETSC_NULL, PETSC_NULL, PETSC_NULL, &max_it );
    // Weirdness with petsc 3.2 here...look at it later
    //BSSCR_KSPLogSetMonitor( ksp_S, max_it, &monitor_index );

    /***************************************************************************************************************/
    /* Pressure / Velocity Solve */     
    /***************************************************************************************************************/
    PetscPrintf( PETSC_COMM_WORLD, "\t* Pressure / Velocity Solve \n");
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    usePreviousGuess = PETSC_FALSE; 
    if(been_here)
        PetscOptionsGetTruth( PETSC_NULL, "-scr_use_previous_guess", &usePreviousGuess, &found );
    if(usePreviousGuess) {   /* Note this should actually look at checkpoint information */
        KSPSetInitialGuessNonzero( ksp_S, PETSC_TRUE ); }
    else {
        KSPSetInitialGuessNonzero( ksp_S, PETSC_FALSE ); }
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /*******     SET CONVERGENCE TESTS     *************************************************************************/
    useNormInfStoppingConditions = PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-scr_use_norm_inf_stopping_condition", &useNormInfStoppingConditions, &found );
    if(useNormInfStoppingConditions) 
        BSSCR_KSPSetNormInfConvergenceTest(ksp_S); 

    useNormInfMonitor = PETSC_FALSE; 
    PetscOptionsGetTruth( PETSC_NULL, "-scr_ksp_norm_inf_monitor", &useNormInfMonitor, &found );
    if(useNormInfMonitor) 
        KSPMonitorSet( ksp_S, BSSCR_KSPNormInfToNorm2Monitor, PETSC_NULL, PETSC_NULL );
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /*******   PRESSURE SOLVE   ************************************************************************************/
    PetscPrintf( PETSC_COMM_WORLD, "\t* KSPSolve( ksp_S, h_hat, p )\n");
    /* if h_hat needs to be fixed up ..take out any nullspace vectors here */
    /* we want to check that there is no "noise" in the null-space in the h vector */
    /* this causes problems when we are trying to solve a Jacobian system when the Residual is almost converged */
    if(bsscrp_self->check_pressureNS){
        bsscrp_self->buildPNS(ksp);/* build and set nullspace vectors on bsscr - which is on ksp (function pointer is set in KSPSetUp_BSSCR) */
    }

    PetscScalar hnorm, gnorm;
    MatNorm(G,NORM_INFINITY,&gnorm);
    VecNorm(h_hat, NORM_2, &hnorm);
    hnorm=hnorm/gnorm;

    if((hnorm < 1e-6) && (hnorm > 1e-20)){
        VecScale(h_hat,1.0/hnorm);
    }
    /* test to see if v or t are in nullspace of G and orthogonalize wrt h_hat if needed */
    KSPRemovePressureNullspace_BSSCR(ksp, h_hat);
    /* set convergence test to use min_it */
    found = PETSC_FALSE;
    min_it = 0;
    PetscOptionsGetInt( PETSC_NULL,"-scr_ksp_set_min_it_converge", &min_it, &found);
    if(found && min_it > 0){
        BSSCR_KSPSetConvergenceMinIts(ksp_S, min_it, bsscrp_self);
    }
    scrSolveTime = MPI_Wtime();
    KSPSolve( ksp_S, h_hat, p );
    scrSolveTime =  MPI_Wtime() - scrSolveTime;
    PetscPrintf( PETSC_COMM_WORLD, "\n\t* KSPSolve( ksp_S, h_hat, p )  Solve  Finished in time: %lf seconds\n\n", scrSolveTime);

    /***************************************/
    if((hnorm < 1e-6) && (hnorm > 1e-20)){
        VecScale(h_hat,hnorm);
        VecScale(p,hnorm);
    }
    KSPRemovePressureNullspace_BSSCR(ksp, p);

    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /* restore K and f for the Velocity back solve */
    found = PETSC_FALSE;
    restorek = PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL, "-restore_K", &restorek, &found);
    //PetscOptionsGetString( PETSC_NULL, "-restore_K", name, PETSC_MAX_PATH_LEN-1, &flg );
    if(penaltyNumber > 1e-10 && bsscrp_self->k2type){
        if(restorek){
            penaltyNumber = -penaltyNumber;
            if(f2) { ierr=VecAXPY(f,penaltyNumber,f2); }/*  f <- f +a*f2 */
            if(f3) { ierr=VecAXPY(f,penaltyNumber,f3); }/*  f <- f +a*f3 */
            ierr=MatAXPY(K,penaltyNumber,K2,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);/* Computes K = penaltyNumber*K2 + K */
            //K=Korig;
            Stg_KSPSetOperators(ksp_inner, K, K, DIFFERENT_NONZERO_PATTERN);
            KisJustK=PETSC_TRUE;          
        }      
    }
    if(f3){ Stg_VecDestroy(&f3 ); } /* always destroy this local vector if was created */

    /* obtain solution for u */
    VecDuplicate( u, &t );   MatMult( G, p, t);  VecAYPX( t, -1.0, f ); /*** t <- -t + f   = f - G*p  ***/
    MatSchurComplementGetKSP( S, &ksp_inner );
    a11SingleSolveTime = MPI_Wtime();           /* ----------------------------------  Final V Solve */
    if(usePreviousGuess) KSPSetInitialGuessNonzero( ksp_inner, PETSC_TRUE );

    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /*******   VELOCITY SOLVE   ************************************************************************************/
    /** Easier to just create a new KSP here if we want to do backsolve diffferently. (getting petsc errors now when switching from fgmres) */
    change_backsolve=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL, "-change_backsolve", &change_backsolve, &found );
    if(change_backsolve){
      //Stg_KSPDestroy(&ksp_inner );
      KSPCreate(PETSC_COMM_WORLD, &ksp_new_inner);
      Stg_KSPSetOperators(ksp_new_inner, K, K, DIFFERENT_NONZERO_PATTERN);
      KSPSetOptionsPrefix(ksp_new_inner, "backsolveA11_");
      KSPSetFromOptions(ksp_new_inner); /* make sure we are setting up our solver how we want it */
      MatSchurComplementSetKSP( S, ksp_new_inner );/* need to give the Schur it's inner ksp back for when we destroy it at end */
      ksp_inner=ksp_new_inner;
    }
    KSPSolve(ksp_inner, t, u);         /* Solve, then restore default tolerance and initial guess */
    a11SingleSolveTime = MPI_Wtime() - a11SingleSolveTime;              /* ------------------ Final V Solve */

    flg=0;
    PetscOptionsGetString( PETSC_NULL, "-solutiondumpdir", name, PETSC_MAX_PATH_LEN-1, &flg );
    if(flg){
        sprintf(str,"%s/",name); sprintf(matname,"p%s",suffix);
        bsscr_dirwriteVec( p, matname,str, "Writing p vector in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"u%s",suffix);
        bsscr_dirwriteVec( u, matname,str, "Writing u vector in al Solver");
        sprintf(str,"%s/",name); sprintf(matname,"h_hat%s",suffix);
        bsscr_dirwriteVec( h_hat, matname,str, "Writing h_hat vector in al Solver");

    }

    found = PETSC_FALSE;
    restorek = PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL, "-restore_K_after_solve", &restorek, &found);
    if(penaltyNumber > 1e-10 && bsscrp_self->k2type){
        if(restorek){
            penaltyNumber = -penaltyNumber;
            ierr=MatAXPY(K,penaltyNumber,K2,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);/* Computes K = penaltyNumber*K2 + K */
            KisJustK=PETSC_TRUE;          
        }      
    }
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    /******     SOLUTION SUMMARY      ******************************************************************************/
    bsscr_summary(bsscrp_self,ksp_S,ksp_inner,K,K2,D,G,C,u,p,f,h,t,penaltyNumber,KisJustK, mgSetupTime, scrSolveTime, a11SingleSolveTime);
    //bsscr_summary(bsscrp_self,ksp_S,ksp_inner,K,Korig,K2,D,G,C,u,p,f,h,t,penaltyNumber,KisJustK, mgSetupTime, scrSolveTime, a11SingleSolveTime);
    /***************************************************************************************************************/
    /***************************************************************************************************************/
    Stg_VecDestroy(&t );
    Stg_KSPDestroy(&ksp_S );
    Stg_VecDestroy(&h_hat );
    Stg_MatDestroy(&S );//This will destroy ksp_inner: also.. pcInner == pc_MG and is destroyed when ksp_inner is 
    been_here = 1;
    PetscFunctionReturn(0);
}
