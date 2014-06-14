#include <mpi.h>
/* petsc.h needs to be included before local .h files for PETSC_VERSION information */
#include <petsc.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>


#include "types.h"
//#include <petscvec.h>
#include "petscvec.h"
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsnes.h>

#include "private/kspimpl.h"   /*I "petscksp.h" I*/

//#include "ksptypes.h"
#include "ksp-register.h"
#include "StokesFullMatrixSNESInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#if ( (PETSC_VERSION_MAJOR>=3) && (PETSC_VERSION_MINOR>=2) )

/* Macro for checking number integrity - i.e. checks if number is infinite or "not a number" */
#define SFMSNES_isGoodNumber( number ) ( (! isnan( number ) ) && ( ! isinf( number ) ) )

#define SFMSNES_GetPetscMatrix( matrix ) ( (Mat)(matrix) )
#define SFMSNES_GetPetscVector( vector ) ( (Vec)(vector) )
#define SFMSNES_GetPetscKSP( solver ) ( (KSP)(solver)  )

        #define __KSP_COMMON \
                Stokes_SLE     *   st_sle; \
                PETScMGSolver  *   mg;	   \
                PetscTruth         DIsSym; \
                StiffnessMatrix*   preconditioner;

        struct KSP_COMMON { __KSP_COMMON };
        typedef struct KSP_COMMON KSP_COMMON;

const Type StokesFullMatrixSNESInterface_Type = "StokesFullMatrixSNESInterface";

void* _StokesFullMatrixSNESInterface_DefaultNew( Name name ) {
    SizeT                                              _sizeOfSelf = sizeof(StokesFullMatrixSNESInterface);
    Type                                                      type = StokesFullMatrixSNESInterface_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLE_Solver_Delete;
    Stg_Class_PrintFunction*                                _print = _SLE_Solver_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLE_Solver_Copy;
    Stg_Component_BuildFunction*                            _build = _StokesFullMatrixSNESInterface_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _StokesFullMatrixSNESInterface_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLE_Solver_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLE_Solver_Destroy;
    SLE_Solver_GetResidualFunc*                       _getResidual = NULL;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StokesFullMatrixSNESInterface_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _StokesFullMatrixSNESInterface_AssignFromXML;
    SLE_Solver_SolverSetupFunction*                   _solverSetup = _StokesFullMatrixSNESInterface_SolverSetup;
    SLE_Solver_SolveFunction*                               _solve = _StokesFullMatrixSNESInterface_Solve;

    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
    return (void*) _StokesFullMatrixSNESInterface_New( STOKESFULLMATRIXSNESINTERFACE_PASSARGS );
}


/* Creation implementation / Virtual constructor */
/* Set up function pointers */
StokesFullMatrixSNESInterface* _StokesFullMatrixSNESInterface_New( STOKESFULLMATRIXSNESINTERFACE_DEFARGS )
{    
    StokesFullMatrixSNESInterface* self;
    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(StokesFullMatrixSNESInterface) );

    self = (StokesFullMatrixSNESInterface*) _SLE_Solver_New( SLE_SOLVER_PASSARGS );

    /* Virtual info */
    return self;
}

void _StokesFullMatrixSNESInterface_Init( 
		StokesFullMatrixSNESInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg )
{
	self->preconditioner = preconditioner;
	self->st_sle = st_sle;
	self->mg     = mg;
}

void _StokesFullMatrixSNESInterface_Build( void* solver, void* sle ) {/* it is the sle here being passed in*/
	StokesFullMatrixSNESInterface*	self  = (StokesFullMatrixSNESInterface*)solver;
	
	Stream_IndentBranch( StgFEM_Debug );
	//if( self->st_sle){
	//  Stg_Component_Build( self->st_sle, dummy, False ); /* causes infinite loop */
	//}

	/* Build Preconditioner */
	if ( self->preconditioner ) {
		Stg_Component_Build( self->preconditioner, sle, False );
		SystemLinearEquations_AddStiffnessMatrix( self->st_sle, self->preconditioner );

	}
	if( self->mg ){
	    Stg_Component_Build( self->mg, sle, False );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}

void _StokesFullMatrixSNESInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data ) {
	StokesFullMatrixSNESInterface* self         = (StokesFullMatrixSNESInterface*) solver;
	//double                  tolerance;
	//Iteration_Index         maxUzawaIterations, minUzawaIterations;
	StiffnessMatrix*        preconditioner;
	//Bool                    useAbsoluteTolerance;
	//Bool                    monitor;
	Stokes_SLE *            st_sle;
	PETScMGSolver *         mg;

	_SLE_Solver_AssignFromXML( self, cf, data );

	preconditioner = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Preconditioner", StiffnessMatrix, False, data  );
	st_sle  = Stg_ComponentFactory_ConstructByName( cf, (Name)"stokesEqn", Stokes_SLE, True, data  ); 
	mg      = Stg_ComponentFactory_ConstructByName( cf, (Name)"mgSolver", PETScMGSolver, False, data);
	self->ctx = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  );

	st_sle->isNonLinear = True;

	_StokesFullMatrixSNESInterface_Init( self, preconditioner, st_sle, mg );

}

void _StokesFullMatrixSNESInterface_Initialise( void* solver, void* stokesSLE ) {
	StokesFullMatrixSNESInterface* self = (StokesFullMatrixSNESInterface*) solver;
	Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
	/* Initialise Parent */
	_SLE_Solver_Initialise( self, sle );
	
	if ( sle->context && (True == sle->context->loadFieldsFromCheckpoint) ) {
		/* The previous timestep's velocity solution will be helpful in iterating to a better
		solution faster - and thus make restarting from checkpoint more repeatable compared
		to original non-restart solution */
		SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->uSolnVec );
		SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->pSolnVec );
	}
	KSPRegisterAllKSP("Solvers/KSPSolvers/src");
}

/* SolverSetup */

void _StokesFullMatrixSNESInterface_SolverSetup( void* solver, void* stokesSLE ) {
	StokesFullMatrixSNESInterface* self = (StokesFullMatrixSNESInterface*) solver;
	//Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Stream_UnIndentBranch( StgFEM_Debug );
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/* to-do: these two functions need to be somewhere common now like a 'Utils' directory */
void sfm_writeMat(Mat K, char name[], char message[]){
    PetscViewer             mat_view_file;
    char str[100];

    if( K != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)K,name);
	sprintf(str,"%smatrixBin",name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_NATIVE );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
	sprintf(str,"%smatrix.m",name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &mat_view_file );
	PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_ASCII_MATLAB );
	MatView( K, mat_view_file );
	Stg_PetscViewerDestroy(&mat_view_file );
    }
}
void sfm_writeVec(Vec V, char name[], char message[]){
    PetscViewer             vec_view_file;
    char str[100];

    if( V != NULL ) {
	PetscPrintf( PETSC_COMM_WORLD,"%s \n",message);
	PetscObjectSetName((PetscObject)V,name);
	sprintf(str,"%svectorBin",name);
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_NATIVE );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
	sprintf(str,"%svector.m",name);
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &vec_view_file );
	PetscViewerSetFormat( vec_view_file, PETSC_VIEWER_ASCII_MATLAB );
	VecView( V, vec_view_file );
	Stg_PetscViewerDestroy(&vec_view_file );
    }
}
/** Assumes A already created and of correct size **/
/** k = (0,1,2,3)
    0 -> insert mat into 0,0 block
    1 -> insert mat into 0,1
    2 -> insert mat into 1,0
    3 -> insert mat into 1,1
    
    K is passed in as it is the stiffness matrix that lives in 0,0 and 
    we can use it to determine the correct index ranges for the other
    matrices being inserted into A.

    (obviously K will equal Asub at some point)
**/
void SFMSNES_InsertSubMatinMat(Mat K, Mat Asub, int k, Mat *A){
    PetscInt m ,n;
    PetscInt km ,kn;
    PetscInt idxm, *idxn;
    MatScalar *a;
    PetscInt          row,ncols;
    const PetscInt    *cols;
    PetscScalar *vals;
    int i,j;
    PetscInt       Istart,Iend;
    PetscInt       rowOffset, colOffset; /* so we can insert matrices into different "blocks" in the larger matrix */

    MatGetSize(K,&km,&kn);
    MatGetSize(Asub,&m,&n);/* m= num rows; n= num cols; */
    rowOffset=(km)*(k/2);
    colOffset=(kn)*(k%2);
    /**    Currently, all PETSc parallel matrix formats are partitioned by
	   contiguous chunks of rows across the processors.  Determine which
	   rows of the matrix are locally owned. **/
    MatGetOwnershipRange(Asub,&Istart,&Iend);
    /* allocate enough memory to hold a full single row of matrix column indices */
    //PetscMalloc(n*sizeof(PetscInt),&idxn);
    for(i=Istart;i<Iend;i++){
	/* ncols = number of non-zero columns for current row */
	/* cols  = indices of non-zero columns */
	/* vals  = array of non-zero values for row */
	MatGetRow(Asub,i,&ncols,&cols,&vals);
	PetscMalloc(ncols*sizeof(PetscInt),&idxn);
	for(j=0;j<ncols;j++){ idxn[j]=cols[j]+colOffset; }
	idxm = i + rowOffset;
	MatSetValues(*A,1,&idxm,ncols,idxn,vals,INSERT_VALUES);
	MatRestoreArray(Asub,&vals);
	PetscFree(idxn);
    }
}
/** Assumes V already created and of correct size **/
/** k = (0,1)
    0 -> insert mat into 0 block 
    1 -> insert mat into 1

    F is passed in as it is the force vector that lives in 0 and 
    we can use it to determine the correct index ranges for the other
    vectors being inserted into V.

**/
void SFMSNES_InsertSubVecinVec(Vec F, Vec Vsub, int k, Vec *V){
    PetscInt n;
    PetscInt kn;
    PetscInt *idxn;
    PetscScalar *a;
    PetscInt          row,ncols;
    //const PetscInt    *cols;
    PetscScalar *vals;
    int i,j;
    PetscInt       Istart,Iend;
    PetscInt       colOffset, ldim;

    VecGetSize(F,&kn);
    VecGetSize(Vsub ,&n);/* n= num rows; */
    VecGetArray(Vsub, &a);

    colOffset=(kn)*(k%2);

    //PetscMalloc(n*sizeof(PetscInt),&idxn);

    VecGetOwnershipRange(Vsub,&Istart,&Iend);
    VecGetLocalSize(Vsub,&ldim);
    PetscMalloc(ldim*sizeof(PetscInt),&idxn);
    PetscMalloc(ldim*sizeof(PetscScalar),&vals);
    for (i=0; i<ldim; i++) {
	idxn[i] = i + Istart;
    }
    VecAssemblyBegin( Vsub );
    VecAssemblyEnd( Vsub );
    VecGetValues(Vsub, ldim, idxn, vals);

    Istart += colOffset;
    for (i=0; i<ldim; i++) {
	idxn[i] += colOffset;;
    }
    VecSetValues(*V,ldim,idxn,vals,INSERT_VALUES);

    VecAssemblyBegin( *V );
    VecAssemblyEnd( *V );

    PetscFree(idxn);
    PetscFree(vals);
}
void SFMSNES_FormMatrixOperator(	Mat A11, Mat A12, Mat A21, Mat A22, Mat* A )
{
    PetscInt km ,kn;
    PetscInt gm, gn;
    if( *A == PETSC_NULL ) {
	MatCreate( MPI_COMM_WORLD, A );
	MatGetSize(A11,&km,&kn);
	MatGetSize(A12,&gm,&gn);
	/* note: kn=km=gm */
	MatSetSizes( *A, PETSC_DECIDE, PETSC_DECIDE, kn+gn, kn+gn );
	MatSetType( *A, MATMPIAIJ );
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
    }

    SFMSNES_InsertSubMatinMat(A11, A11, 0, A);
    SFMSNES_InsertSubMatinMat(A11, A12, 1, A);
    SFMSNES_InsertSubMatinMat(A11, A21, 2, A);
    if( A22 ) SFMSNES_InsertSubMatinMat(A11, A22, 3, A);

    /* if( A11 ) MatBlockSetValue( *A, 0, 0, A11, SAME_NONZERO_PATTERN, INSERT_VALUES ); */
    /* if( A12 ) MatBlockSetValue( *A, 0, 1, A12, SAME_NONZERO_PATTERN, INSERT_VALUES ); */
    /* if( A21 ) MatBlockSetValue( *A, 1, 0, A21, SAME_NONZERO_PATTERN, INSERT_VALUES ); */
    /* if( A22 ) MatBlockSetValue( *A, 1, 1, A22, SAME_NONZERO_PATTERN, INSERT_VALUES ); */

    MatAssemblyBegin( *A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( *A, MAT_FINAL_ASSEMBLY );
}
/* returns A x and b regular (non-block) Petsc Mat and Vecs */
void SFMSNES_FormMatrixSystem( Mat A11, Mat A12, Mat A21, Mat A22,
			     Vec u, Vec p, Vec f, Vec h,
			     Mat* A, Vec *x, Vec *b )
{
   Vec uDup, pDup, fDup, hDup;

   SFMSNES_FormMatrixOperator( A11, A12, A21, A22, A );

   MatGetVecs( *A, x, b );
   VecDuplicate( u, &uDup ); VecCopy( u, uDup );
   VecDuplicate( p, &pDup ); VecCopy( p, pDup );
   VecDuplicate( f, &fDup ); VecCopy( f, fDup );
   VecDuplicate( h, &hDup ); VecCopy( h, hDup );

   SFMSNES_InsertSubVecinVec(uDup, uDup, 0, x);
   SFMSNES_InsertSubVecinVec(uDup, pDup, 1, x);

   SFMSNES_InsertSubVecinVec(fDup, fDup, 0, b);
   SFMSNES_InsertSubVecinVec(fDup, hDup, 1, b);

   Stg_VecDestroy(&uDup );
   Stg_VecDestroy(&pDup );
   Stg_VecDestroy(&fDup );
   Stg_VecDestroy(&hDup );

   VecAssemblyBegin( *x );
   VecAssemblyBegin( *b );
   VecAssemblyEnd( *x );
   VecAssemblyEnd( *b );
}

/* get Mats and Vecs from the SLE */
void SFMSNES_GetStokesOperators( 
		Stokes_SLE *stokesSLE,
		Mat *K,Mat *G,Mat *D,Mat *C,Mat *Smat,
		Vec *f,Vec *h,Vec *u,Vec *p )
{
	
    if(K) {if (stokesSLE->kStiffMat){      *K = SFMSNES_GetPetscMatrix( stokesSLE->kStiffMat->matrix );     } else { *K = PETSC_NULL;}}
    if(G) {if (stokesSLE->gStiffMat){      *G = SFMSNES_GetPetscMatrix( stokesSLE->gStiffMat->matrix );     } else { *G = PETSC_NULL;}}
    if(D) {if (stokesSLE->dStiffMat){      *D = SFMSNES_GetPetscMatrix( stokesSLE->dStiffMat->matrix );     } else { *D = PETSC_NULL;}}
    if(C) {if (stokesSLE->cStiffMat){      *C = SFMSNES_GetPetscMatrix( stokesSLE->cStiffMat->matrix );     } else { *C = PETSC_NULL;}}
	
	/* preconditioner */
	if(Smat){
	    *Smat = PETSC_NULL;
	    if( ((StokesFullMatrixSNESInterface*)stokesSLE->solver)->preconditioner ) {
		StiffnessMatrix *preconditioner;
		
		preconditioner = ((StokesFullMatrixSNESInterface*)stokesSLE->solver)->preconditioner;
		*Smat = SFMSNES_GetPetscMatrix( preconditioner->matrix );
	    }
	}
	if(f){
	    *f = PETSC_NULL;
	    if (stokesSLE->fForceVec){      *f = SFMSNES_GetPetscVector( stokesSLE->fForceVec->vector );     }
	}
	if(h){
	    *h = PETSC_NULL;
	    if (stokesSLE->hForceVec){      *h = SFMSNES_GetPetscVector( stokesSLE->hForceVec->vector );     }
	}
	if(u){
	    *u = PETSC_NULL;
	    if (stokesSLE->uSolnVec){       *u = SFMSNES_GetPetscVector( stokesSLE->uSolnVec->vector );      }
	}
	if(p){
	    *p = PETSC_NULL;
	    if (stokesSLE->pSolnVec){       *p = SFMSNES_GetPetscVector( stokesSLE->pSolnVec->vector );      }
	}
}
void SFMSNES_Assemble( StokesFullMatrixSNESInterface* self, Bool nonLinear ) {
   int index;
   PetscTruth printMats, found;
   //SystemLinearEquations* SLE = (SystemLinearEquations*)(self->st_sle);
   Stokes_SLE* SLE = (Stokes_SLE*)(self->st_sle);
   
   SystemLinearEquations_UpdateSolutionOntoNodes( self->st_sle, self->ctx );
   SystemLinearEquations_ZeroAllVectors( self->st_sle, self->ctx );/* zeroes all force vectors only */
   printMats = PETSC_FALSE;
   PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_printmat", &printMats, &found );

   self->st_sle->nlFormJacobian = nonLinear;
   for ( index = 0; index < SLE->forceVectors->count; index++ ) {
       ForceVector_Assemble( SLE->forceVectors->data[index] );
       if(printMats){
	   PetscPrintf( PETSC_COMM_WORLD, "Force vector name %s\n",SLE->forceVectors->data[index]->name); }
   }
   for ( index = 0; index < SLE->stiffnessMatrices->count; index++ ) {
       StiffnessMatrix_Assemble( SLE->stiffnessMatrices->data[index], SLE->removeBCs, SLE, self->ctx );
       if(printMats){
	   PetscPrintf( PETSC_COMM_WORLD, "Stiffness matrix name %s\n",SLE->stiffnessMatrices->data[index]->name); }
   }
   /* StiffnessMatrix_Assemble( SLE->kStiffMat, SLE->removeBCs, SLE, self->ctx );  */
   /* StiffnessMatrix_Assemble( SLE->gStiffMat, SLE->removeBCs, SLE, self->ctx ); */
   /* PetscPrintf( PETSC_COMM_WORLD, "K and G Stokes matrices re-assembled in %s\n",__func__); */
   if(printMats){
       PetscPrintf( PETSC_COMM_WORLD, ":: num Matrices = %d :: num Force Vec = %d \n", SLE->stiffnessMatrices->count, SLE->forceVectors->count); }
   self->st_sle->nlFormJacobian = False;
}
PetscErrorCode SFMSNES_Check(Mat M, Vec V, const char *f, char* Mname, char* Vname){
    PetscReal norm;
    int m,n;
    PetscTruth assembled;
    if(M){
    	MatGetSize( M, &m, &n);
    	MatNorm(M,NORM_1,&norm);
    	PetscPrintf( PETSC_COMM_WORLD,  "\t\t Mat %s norm is %.6e from %s: address is %p: size : %dx%d\n",Mname, norm, f ,M, m,n);
    	MatAssembled( M, &assembled);
    	if(assembled){ 	PetscPrintf( PETSC_COMM_WORLD,  "\t\t Mat is assembled\n"); }
    	else{ PetscPrintf( PETSC_COMM_WORLD,  "\t\t Mat is NOT assembled\n"); }
    }
    if(V){
    	VecNorm( V, NORM_2, &norm );
    	VecGetSize(V,&n);
    	PetscPrintf( PETSC_COMM_WORLD,  "\t\t Vec %s norm is %.6e from %s: address is %p: size : %d\n", Vname, norm, f ,V, n);
    	//if(V->assembled){ 	PetscPrintf( PETSC_COMM_WORLD,  "\t\t Vec is assembled\n"); }
    	//else{ PetscPrintf( PETSC_COMM_WORLD,  "\t\t Vec is NOT assembled\n"); }
    }
    return 0;
}
#include <private/snesimpl.h>    /*I  "petscsnes.h"  I*/

#undef __FUNCT__  
#define __FUNCT__ "SFMSSNES_FDComputeJacobian"
PetscErrorCode  SFMSSNES_FDComputeJacobian(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
  Vec            j1a,j2a,x2;
  PetscErrorCode ierr;
  PetscInt       i,N,start,end,j,value,root;
  PetscScalar    dx,*y,*xx,wscale;
  PetscReal      amax,epsilon = PETSC_SQRT_MACHINE_EPSILON;
  PetscReal      dx_min = 1.e-16,dx_par = 1.e-1,unorm;
  MPI_Comm       comm;
  PetscErrorCode (*eval_fct)(SNES,Vec,Vec)=0;
  PetscBool      assembled,use_wp = PETSC_TRUE,flg;
  const char     *list[2] = {"ds","wp"};
  PetscMPIInt    size;
  const PetscInt *ranges;

  PetscFunctionBegin;
  ierr = PetscOptionsGetReal(((PetscObject)snes)->prefix,"-snes_test_err",&epsilon,0);CHKERRQ(ierr);
  eval_fct = SNESComputeFunction;

  ierr = PetscObjectGetComm((PetscObject)x1,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MatAssembled(*B,&assembled);CHKERRQ(ierr);
  if (assembled) {
    ierr = MatZeroEntries(*B);CHKERRQ(ierr);
  }
  if (!snes->nvwork) {
    snes->nvwork = 3;
    ierr = VecDuplicateVecs(x1,snes->nvwork,&snes->vwork);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(snes,snes->nvwork,snes->vwork);CHKERRQ(ierr);
  }
  j1a = snes->vwork[0]; j2a = snes->vwork[1]; x2 = snes->vwork[2];

  ierr = VecGetSize(x1,&N);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(x1,&start,&end);CHKERRQ(ierr);
  ierr = (*eval_fct)(snes,x1,j1a);CHKERRQ(ierr);

  ierr = PetscOptionsEList("-mat_fd_type","Algorithm to compute difference parameter","SNESDefaultComputeJacobian",list,2,"wp",&value,&flg);CHKERRQ(ierr);
  if (flg && !value) {
    use_wp = PETSC_FALSE;
  }
  if (use_wp) {
    ierr = VecNorm(x1,NORM_2,&unorm);CHKERRQ(ierr);
  }
  /* Compute Jacobian approximation, 1 column at a time. 
      x1 = current iterate, j1a = F(x1)
      x2 = perturbed iterate, j2a = F(x2)
   */
  for (i=0; i<N; i++) {
    ierr = VecCopy(x1,x2);CHKERRQ(ierr);
    if (i>= start && i<end) {
      ierr = VecGetArray(x1,&xx);CHKERRQ(ierr);
      if (use_wp) {
        dx = 1.0 + unorm;
      } else {
        dx = xx[i-start];
      }
      ierr = VecRestoreArray(x1,&xx);CHKERRQ(ierr);
      if (PetscAbsScalar(dx) < dx_min) dx = (PetscRealPart(dx) < 0. ? -1. : 1.) * dx_par;
      dx *= epsilon;
      wscale = 1.0/dx;
      ierr = VecSetValues(x2,1,&i,&dx,ADD_VALUES);CHKERRQ(ierr);
    } else {
      wscale = 0.0;
    }
    ierr = VecAssemblyBegin(x2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x2);CHKERRQ(ierr);
    ierr = (*eval_fct)(snes,x2,j2a);CHKERRQ(ierr);
    ierr = VecAXPY(j2a,-1.0,j1a);CHKERRQ(ierr);
    /* Communicate scale=1/dx_i to all processors */
    ierr = VecGetOwnershipRanges(x1,&ranges);CHKERRQ(ierr);
    root = size;
    for (j=size-1; j>-1; j--){
      root--;
      if (i>=ranges[j]) break;
    }
    ierr = MPI_Bcast(&wscale,1,MPIU_SCALAR,root,comm);CHKERRQ(ierr);

    ierr = VecScale(j2a,wscale);CHKERRQ(ierr);
    ierr = VecNorm(j2a,NORM_INFINITY,&amax);CHKERRQ(ierr); amax *= 1.e-14;
    ierr = VecGetArray(j2a,&y);CHKERRQ(ierr);
    for (j=start; j<end; j++) {
      if (PetscAbsScalar(y[j-start]) > amax || j == i) {
        ierr = MatSetValues(*B,1,&j,1,&i,y+j-start,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecRestoreArray(j2a,&y);CHKERRQ(ierr);
  }
  ierr  = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr  = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*B != *J) {
    ierr  = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr  = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  *flag =  DIFFERENT_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}
/* w is a normal Petsc vector that contains current guess from Newton SNES method */
/* to get Residual we need to get the u and p "parts" of it and update u and p onto mesh
   and re-assemble A(w)*w-b(w) = F(w) */
/* this function stomps on the f and h force vectors on the SLE as well */
PetscErrorCode SFMSNES_FormResidual( SNES snes, Vec w, Vec F, StokesFullMatrixSNESInterface* self )
{
   Mat K, G, D, C, S;
   Vec u, p, f, h, *subVecs;
   Vec save_u, save_p;/* temporary vectors to hold current solution */
   Vec t,v; /* temp vecs to hold null-space vectors */
   double wallTime, uNorm, pNorm, bNorm;
   static double totTime=0.0;
   static int timesCalled=0;
   PetscTruth useTimer, flg, found, uDiff=PETSC_FALSE;
   Stokes_SLE* SLE = (Stokes_SLE*)(self->st_sle);
   IJK          elementRes;
   Index nx,ny,nz;
   PetscReal eps = PETSC_SQRT_MACHINE_EPSILON;

   IS is;
   PetscInt un, pn;
   PetscInt j,start,end;
   int i;
   Vec subV;
   PetscInt       *indices;
   PetscErrorCode          ierr;

   useTimer = PETSC_FALSE;
   PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_use_timer", &useTimer, &found );
   if(useTimer){
       wallTime = MPI_Wtime();
   }

   /* w is current solution being iterated on */
   /* u and p are the SLE solution vectors */
   /* so here we are setting the SLE vectors to the current solution */
   /* update solution onto nodes  happens in the assemble function */
   SFMSNES_GetStokesOperators( self->st_sle, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL,  PETSC_NULL, PETSC_NULL, &u, &p );
   /* SFMSNES_Check(0,u, __func__, 0, "u"); */
   /* SFMSNES_Check(0,p, __func__, 0, "p"); */

   /* we want to overwrite u and p with the values in w */
   /* w is a "full" size vector w=[u_i;p_i] i is ith iteration */
   VecGetSize(u ,&un);
   PetscMalloc(un*sizeof(PetscInt),&indices);
   for(i=0;i<un;i++){
       indices[i]=i;
   }
   ISCreateGeneral(PETSC_COMM_SELF,un,indices, PETSC_COPY_VALUES,&is);
   VecGetSubVector( w, is, &subV);
   VecCopy(subV, u);/* stomp on u vector now */
   VecRestoreSubVector( w, is, &subV);
   PetscFree(indices);
   Stg_ISDestroy(&is);

   VecGetSize(p ,&pn);
   PetscMalloc(pn*sizeof(PetscInt),&indices);
   for(i=0;i<pn;i++){
       indices[i]=i+un;
   }
   ISCreateGeneral(PETSC_COMM_SELF,pn,indices, PETSC_COPY_VALUES,&is);
   VecGetSubVector( w, is, &subV);
   VecCopy(subV, p);/* stomp on p vector now */
   VecRestoreSubVector( w, is, &subV);
   PetscFree(indices);
   Stg_ISDestroy(&is);

   /* update Mats and Vecs with new values of u and p */
   SFMSNES_Assemble( self, False );

   SFMSNES_GetStokesOperators( self->st_sle, &K, &G, &D, &C, PETSC_NULL, &f, &h, PETSC_NULL, PETSC_NULL );

   VecScale( f, -1.0 );
   VecScale( h, -1.0 );
   MatMultAdd( K, u, f, f ); /* K*u + f -> f */
   MatMultAdd( G, p, f, f );
   if( D ) MatMultAdd( D, u, h, h );
   else MatMultTransposeAdd( G, u, h, h );
   if( C ) MatMultAdd( C, p, h, h );

   /* we want to check that there is no "noise" in the null-space in the h vector */
   /* this causes problems when we are trying to solve a Jacobian system when the Residual is almost converged */
   /* first construct checkerboard part of possible null-space */
   /* Get Mesh Sizes first */
   nx=elementRes[I_AXIS] = Dictionary_GetInt( self->ctx->CF->rootDict, (Dictionary_Entry_Key)"elementResI"  );
   ny=elementRes[J_AXIS] = Dictionary_GetInt( self->ctx->CF->rootDict, (Dictionary_Entry_Key)"elementResJ"  );
   nz=elementRes[K_AXIS] = Dictionary_GetInt( self->ctx->CF->rootDict, (Dictionary_Entry_Key)"elementResK"  );

   VecDuplicate(h,&t);
   VecGetOwnershipRange(t,&start,&end);
   /* for 2D atm */
   for (j=start; j<end; j++) {
        if( j%2 == (j/nx)%2 ){
	    VecSetValue(t,j,1.0,INSERT_VALUES);
	}
   }
   VecAssemblyBegin( t );
   VecAssemblyEnd( t );

   VecDuplicate(h,&v);
   VecGetOwnershipRange(v,&start,&end);
   /* for 2D atm */
   for (j=start; j<end; j++) {
        if( j%2 != (j/nx)%2 ){
	    VecSetValue(v,j,1.0,INSERT_VALUES);
	}
   }
   VecAssemblyBegin( v );
   VecAssemblyEnd( v );
   PetscScalar norm, a, a1, a2;
   Vec t2;
   eps=1e-19;
   /* test to see if v or t are in nullspace of G */
   VecDuplicate(u,&t2);
   MatMult( G, t, t2);
   VecNorm(t2, NORM_2, &norm);
   if(norm < 1e-10){/* then t in NS of G */
       VecNorm(h, NORM_2, &norm);
       if(norm > 10*eps){ VecScale(h,1.0/norm); }
       VecDot(t,h, &a1);
       VecDot(t,t, &a2);
       a=-a1/a2;
       VecAXPY(h, a, t);
       VecDot(t,h, &a1);
       if(norm > 10*eps){ VecScale(h,norm); }
   }
   MatMult( G, v, t2);
   VecNorm(t2, NORM_2, &norm);
   if(norm < 1e-10){/* then t in NS of G */
       VecNorm(h, NORM_2, &norm);
       if(norm > 10*eps){VecScale(h,1.0/norm); }
       VecDot(v,h, &a1);
       VecDot(v,v, &a2);
       a=-a1/a2;
       VecAXPY(h, a, v);
       VecDot(v,h, &a1);
       if(norm > 10*eps){ VecScale(h,norm); }
   }
   VecDot(t,h, &a1);
   VecDot(v,h, &a1);
/*
  a=(v'*r)/(v'*v);
  b=(t'*r)/(t'*t);

  rr=r-a*v-b*t;
*/
/*
  if (mod(i-1,2) == 0) && (mod(floor((i-1)/nEx),2)==0 )
    t(i)=1;
  endif
  if (mod(i-1,2) == 1) && (mod(floor((i-1)/nEx),2)==1 )
    t(i)=1;
  endif
*/

   SFMSNES_InsertSubVecinVec(f, f, 0, &F);
   SFMSNES_InsertSubVecinVec(f, h, 1, &F);
   VecAssemblyBegin( F );
   VecAssemblyEnd( F );

   flg = PETSC_FALSE;
   ierr = PetscOptionsGetBool(PETSC_NULL,"-snes_dump_residual",&flg,PETSC_NULL);CHKERRQ(ierr);
   if(flg){
       sfm_writeVec( F, "F", "Dumping Residual Vector");
       sfm_writeVec( t, "t", "Dumping t NS Vector");
       sfm_writeVec( v, "v", "Dumping v NS Vector");
       sfm_writeVec( h, "H", "Dumping H Vector");
   }

   if(useTimer){
       wallTime = MPI_Wtime() - wallTime;
       PetscPrintf( PETSC_COMM_WORLD, "  F assemble time = %f \n", wallTime );
       totTime += wallTime;
       timesCalled++;
       PetscPrintf( PETSC_COMM_WORLD, "  F assemble time = %f : Accumulated Total. Called %d times \n", totTime, timesCalled );
   }

   return 0;
}
/* J and Jp get passed back to SNES and are used as the mat and pmat on the PC in the KSP */
PetscErrorCode SFMSNES_FormJacobian(
   SNES snes, Vec x, Mat* J, Mat* Jp, MatStructure* flag,
   StokesFullMatrixSNESInterface* self )
{
    Mat K, G, D, C, Gt, Shat;
    Vec u, p, *subVecs;
    Vec save_u, save_p;/* temporary vectors to hold current solution */
    double wallTime;
    static double totTime=0.0;
    static int timesCalled=0;
    PetscTruth useTimer, flg, found, uDiff=PETSC_FALSE;
    
    IS is;
    PetscInt un, pn;
    int i;
    Vec subV;
    PetscInt       *indices;
    PetscErrorCode          ierr;
    MatStructure mstruct;

    useTimer = PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_use_timer", &useTimer, &found );
    if(useTimer){
	wallTime = MPI_Wtime();
    }

    SFMSNES_GetStokesOperators( self->st_sle, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL,
			     &u, &p );
    /* x is current solution being iterated on */
    /* u and p are the SLE solution vectors */
    /* SFMSNES_Check(0,u, __func__, 0, "u"); */
    /* SFMSNES_Check(0,p, __func__, 0, "p"); */
    /* we want to overwrite u and p with the values in w */
   /* w is a "full" size vector w=[u_i;p_i] i is ith iteration */
   VecGetSize(u ,&un);
   PetscMalloc(un*sizeof(PetscInt),&indices);
   for(i=0;i<un;i++){
       indices[i]=i;
   }
   ISCreateGeneral(PETSC_COMM_SELF,un,indices, PETSC_COPY_VALUES,&is);
   VecGetSubVector( x, is, &subV);
   VecCopy(subV, u);/* stomp on u vector now */
   VecRestoreSubVector( x, is, &subV);
   PetscFree(indices);
   Stg_ISDestroy(&is);

   VecGetSize(p ,&pn);
   PetscMalloc(pn*sizeof(PetscInt),&indices);
   for(i=0;i<pn;i++){
       indices[i]=i+un;
   }
   ISCreateGeneral(PETSC_COMM_SELF,pn,indices, PETSC_COPY_VALUES,&is);
   VecGetSubVector( x, is, &subV);
   VecCopy(subV, p);/* stomp on p vector now */
   VecRestoreSubVector( x, is, &subV);
   PetscFree(indices);
   Stg_ISDestroy(&is);

   /* update Mats and Vecs with new values of u and p */
   SFMSNES_Assemble( self, True );
   SFMSNES_GetStokesOperators( self->st_sle, &K, &G, &D, &C, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL );
   /* create Gt if no D */
   if( !D ) {
       MatTranspose(G, MAT_INITIAL_MATRIX,&Gt);
   }
   else {
       Gt = D;
   }

   flg = PETSC_FALSE;
   ierr = PetscOptionsGetBool(PETSC_NULL,"-snes_dump_suboperators",&flg,PETSC_NULL);CHKERRQ(ierr);
   if(flg){
       sfm_writeMat( K, "K", "Dumping K Matrix");
       sfm_writeMat( G, "G", "Dumping G Matrix");
       sfm_writeMat( Gt, "Gt", "Dumping Gt Matrix");
       sfm_writeMat( C, "C", "Dumping C Matrix");
   }

   SFMSNES_FormMatrixOperator( K,G,Gt,C, J );
   flg = PETSC_FALSE;
   ierr = PetscOptionsGetBool(PETSC_NULL,"-snes_dump_operator",&flg,PETSC_NULL);CHKERRQ(ierr);
   if(flg){
       sfm_writeMat( *J, "J", "Dumping Jacobian Matrix");
   }

   if( !D ) { Stg_MatDestroy(&Gt ); }

   flg = PETSC_FALSE;
   ierr = PetscOptionsGetBool(PETSC_NULL,"-snes_dump_fdJ",&flg,PETSC_NULL);CHKERRQ(ierr);
   if(flg){
       Mat FDexp;

       ierr = MatConvert(*J,MATSAME,MAT_INITIAL_MATRIX,&FDexp);CHKERRQ(ierr);
       ierr = SFMSSNES_FDComputeJacobian(snes,x,&FDexp,&FDexp,&mstruct,NULL);CHKERRQ(ierr);
       
       sfm_writeMat( FDexp, "Jfd", "Dumping FDJacobian Matrix");
   }

   if(useTimer){
       wallTime = MPI_Wtime() - wallTime;
       PetscPrintf( PETSC_COMM_WORLD, "  J assemble time = %f \n", wallTime );
       totTime += wallTime;
       timesCalled++;
       PetscPrintf( PETSC_COMM_WORLD, "  J assemble time = %f : Accumulated Total. Called %d times \n", totTime, timesCalled );
   }

   return 0;
}

#undef __FUNCT__  
#define __FUNCT__ "SFMSNES_DumpMat"
PetscErrorCode SFMSNES_DumpMat(KSP ksp,PetscInt n,PetscReal rnorm, void *dummy)
{
      PetscErrorCode          ierr;
      PetscViewer mat_view_file;
      Mat         A,B;
      char name[10];
      char str[20];
      
      PetscFunctionBegin;
      
      sprintf(name,"%s","A");
      ierr = PCGetOperators(ksp->pc,&A,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
      ierr = MatComputeExplicitOperator(A,&B);CHKERRQ(ierr);
      PetscObjectSetName((PetscObject)B,name);

      sprintf(str,"%smatrixBin",name);
      PetscViewerBinaryOpen( PETSC_COMM_WORLD, str, FILE_MODE_WRITE, &mat_view_file );
      PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_NATIVE );
      MatView( B, mat_view_file );
      PetscViewerDestroy(&mat_view_file );

      sprintf(str,"%smatrix.m",name);
      PetscViewerASCIIOpen( PETSC_COMM_WORLD, str, &mat_view_file );
      PetscViewerSetFormat( mat_view_file, PETSC_VIEWER_ASCII_MATLAB );
      MatView( B, mat_view_file );
      PetscViewerDestroy(&mat_view_file );

      ierr = MatDestroy(&B);CHKERRQ(ierr);
	
      PetscFunctionReturn(0);
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/* Sets up Solver to be a SNES with Petsc ksp (gmres) solve by default: */
void _StokesFullMatrixSNESInterface_Solve( void* solver, void* _stokesSLE ) {
        Stokes_SLE*  stokesSLE  = (Stokes_SLE*)_stokesSLE;
	StokesFullMatrixSNESInterface* self    = (StokesFullMatrixSNESInterface*)solver;
	PetscErrorCode ierr;
	/* Create shortcuts to stuff needed on sle */
	Mat       K;
	Mat       G;
	Mat       Gt;
	Mat       D;
	Mat       C;
	Mat       Smat;
	Vec       u;
	Vec       p;
	Vec       f;
	Vec       h;
	Vec       F;

	Stream*   errorStream = Journal_Register( Error_Type, (Name)StokesFullMatrixSNESInterface_Type  );

	Mat stokes_P;
	Mat stokes_A=0;
	Vec stokes_x;
	Vec stokes_b;

	SNES snes;
	KSP stokes_ksp;
	PC  stokes_pc;
	char name[100];
	
	PetscTruth sym, flg;


	SFMSNES_GetStokesOperators( stokesSLE, &K,&G,&D,&C, &Smat, &f,&h, &u,&p );
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetBool(PETSC_NULL,"-snes_dump_suboperators2",&flg,PETSC_NULL);CHKERRQ(ierr);
	if(flg){
	    sfm_writeMat( K, "K", "Dumping K Matrix");
	    sfm_writeMat( G, "G", "Dumping G Matrix");
	}

	/* create Gt if no D */
	if( !D ) {
	    MatTranspose(G, MAT_INITIAL_MATRIX,&Gt);
	}
	else {
	    Gt = D;
	}
	sym = PETSC_FALSE;

	SFMSNES_FormMatrixSystem( K,G,Gt,C, u,p, f,h, &stokes_A, &stokes_x, &F );

	SNESCreate( PETSC_COMM_WORLD, &snes );

	SNESSetJacobian( snes, stokes_A, stokes_A, SFMSNES_FormJacobian, self );
	SNESSetFunction( snes, F, SFMSNES_FormResidual, self );

	SNESGetKSP( snes, &stokes_ksp );
	KSPSetType( stokes_ksp, "gmres" );/* making this the default solver */
	KSPGetPC( stokes_ksp, &stokes_pc );
	PCSetType( stokes_pc, PCNONE );
	KSPSetFromOptions( stokes_ksp );

	/* flg = PETSC_FALSE; */
	/* ierr = PetscOptionsGetBool(((PetscObject)stokes_ksp)->prefix,"-ksp_dump_operator",&flg,PETSC_NULL);CHKERRQ(ierr); */
	/* if(flg)  KSPMonitorSet( stokes_ksp, SFMSNES_DumpMat, PETSC_NULL, PETSC_NULL ); */

	SNESSetFromOptions( snes );
	SNESSolve( snes, PETSC_NULL, stokes_x );

	Stg_SNESDestroy(&snes );

	Stg_MatDestroy(&stokes_A );
	Stg_VecDestroy(&stokes_x );
	Stg_VecDestroy(&F );

	if(!D){ Stg_MatDestroy(&Gt); }
}

#endif
