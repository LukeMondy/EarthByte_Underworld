#if 0
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsnes.h>
#include <petscsys.h>

#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include "petsc-private/kspimpl.h"   /*I "petscksp.h" I*/
#else
  #include "private/kspimpl.h"   /*I "petscksp.h" I*/
#endif

//#include "ksptypes.h"
#include "ksp-register.h"
#include "StokesBlockSNESInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Macro for checking number integrity - i.e. checks if number is infinite or "not a number" */
#define SBSNES_isGoodNumber( number ) ( (! isnan( number ) ) && ( ! isinf( number ) ) )

#define SBSNES_GetPetscMatrix( matrix ) ( (Mat)(matrix) )
#define SBSNES_GetPetscVector( vector ) ( (Vec)(vector) )
#define SBSNES_GetPetscKSP( solver ) ( (KSP)(solver)  )

        #define __KSP_COMMON \
                Stokes_SLE     *   st_sle; \
                PETScMGSolver  *   mg;	   \
                PetscTruth         DIsSym; \
                StiffnessMatrix*   preconditioner;

        struct KSP_COMMON { __KSP_COMMON };
        typedef struct KSP_COMMON KSP_COMMON;

const Type StokesBlockSNESInterface_Type = "StokesBlockSNESInterface";

void* _StokesBlockSNESInterface_DefaultNew( Name name ) {
    SizeT                                              _sizeOfSelf = sizeof(StokesBlockSNESInterface);
    Type                                                      type = StokesBlockSNESInterface_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLE_Solver_Delete;
    Stg_Class_PrintFunction*                                _print = _SLE_Solver_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLE_Solver_Copy;
    Stg_Component_BuildFunction*                            _build = _StokesBlockSNESInterface_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _StokesBlockSNESInterface_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLE_Solver_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLE_Solver_Destroy;
    SLE_Solver_GetResidualFunc*                       _getResidual = NULL;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StokesBlockSNESInterface_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _StokesBlockSNESInterface_AssignFromXML;
    SLE_Solver_SolverSetupFunction*                   _solverSetup = _StokesBlockSNESInterface_SolverSetup;
    SLE_Solver_SolveFunction*                               _solve = _StokesBlockSNESInterface_Solve;

    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
    return (void*) _StokesBlockSNESInterface_New( STOKESBLOCKSNESINTERFACE_PASSARGS );
}


/* Creation implementation / Virtual constructor */
/* Set up function pointers */
StokesBlockSNESInterface* _StokesBlockSNESInterface_New( STOKESBLOCKSNESINTERFACE_DEFARGS )
{    
    StokesBlockSNESInterface* self;
    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(StokesBlockSNESInterface) );

    self = (StokesBlockSNESInterface*) _SLE_Solver_New( SLE_SOLVER_PASSARGS );

    /* Virtual info */
    return self;
}

void _StokesBlockSNESInterface_Init( 
		StokesBlockSNESInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg )
{
	self->preconditioner = preconditioner;
	self->st_sle = st_sle;
	self->mg     = mg;
}

void _StokesBlockSNESInterface_Build( void* solver, void* sle ) {/* it is the sle here being passed in*/
	StokesBlockSNESInterface*	self  = (StokesBlockSNESInterface*)solver;
	
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

void _StokesBlockSNESInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data ) {
	StokesBlockSNESInterface* self         = (StokesBlockSNESInterface*) solver;
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

	_StokesBlockSNESInterface_Init( self, preconditioner, st_sle, mg );

}

void _StokesBlockSNESInterface_Initialise( void* solver, void* stokesSLE ) {
	StokesBlockSNESInterface* self = (StokesBlockSNESInterface*) solver;
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

void _StokesBlockSNESInterface_SolverSetup( void* solver, void* stokesSLE ) {
	StokesBlockSNESInterface* self = (StokesBlockSNESInterface*) solver;
	//Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Stream_UnIndentBranch( StgFEM_Debug );
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
void SBSNES_FormBlockOperator(	Mat A11, Mat A12, Mat A21, Mat A22, Mat* A )
{
    if( *A == PETSC_NULL ) {
	MatCreate( MPI_COMM_WORLD, A );
	MatSetSizes( *A, 2, 2, 2, 2 );
	MatSetType( *A, "block" );
#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=3)) || (PETSC_VERSION_MAJOR>3) )
        MatSetUp(*A);
#endif
#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=3)) || (PETSC_VERSION_MAJOR>3) )
	MatSetSizes_Block( *A, 2, 2, 2, 2 );
#endif
    }

    if( A11 ) MatBlockSetValue( *A, 0, 0, A11, SAME_NONZERO_PATTERN, INSERT_VALUES );
    if( A12 ) MatBlockSetValue( *A, 0, 1, A12, SAME_NONZERO_PATTERN, INSERT_VALUES );
    if( A21 ) MatBlockSetValue( *A, 1, 0, A21, SAME_NONZERO_PATTERN, INSERT_VALUES );
    if( A22 ) MatBlockSetValue( *A, 1, 1, A22, SAME_NONZERO_PATTERN, INSERT_VALUES );

    MatAssemblyBegin( *A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( *A, MAT_FINAL_ASSEMBLY );
}
void SBSNES_FormBlockSystem( Mat A11, Mat A12, Mat A21, Mat A22,
			     Vec u, Vec p, Vec f, Vec h,
			     Mat* A, Vec *x, Vec *b )
{
   Vec uDup, pDup, fDup, hDup;

   SBSNES_FormBlockOperator( A11, A12, A21, A22, A );

   MatGetVecs( *A, x, b );
   VecDuplicate( u, &uDup ); VecCopy( u, uDup );
   VecDuplicate( p, &pDup ); VecCopy( p, pDup );
   VecDuplicate( f, &fDup ); VecCopy( f, fDup );
   VecDuplicate( h, &hDup ); VecCopy( h, hDup );
   VecBlockSetValue( *x, 0, uDup, INSERT_VALUES );
   VecBlockSetValue( *x, 1, pDup, INSERT_VALUES );
   VecBlockSetValue( *b, 0, fDup, INSERT_VALUES );
   VecBlockSetValue( *b, 1, hDup, INSERT_VALUES );
   Stg_VecDestroy(&uDup );
   Stg_VecDestroy(&pDup );
   Stg_VecDestroy(&fDup );
   Stg_VecDestroy(&hDup );

   VecAssemblyBegin( *x );
   VecAssemblyBegin( *b );
   VecAssemblyEnd( *x );
   VecAssemblyEnd( *b );
}

/** currently not being used **/
PetscErrorCode SBSNES_CreateStokesBlockOperators( MPI_Comm comm, 
					   Mat K, Mat G, Mat D, Mat C,
					   Vec u, Vec p, Vec f, Vec h,
					   Mat *A, Vec *x, Vec *b )
{
    
    MatCreate( comm, A );
    MatSetType( *A, "block" );
    MatSetSizes( *A, 2,2, 2,2 );
#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=3)) || (PETSC_VERSION_MAJOR>3) )
    MatSetUp(*A);
#endif
#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=3)) || (PETSC_VERSION_MAJOR>3) )
    MatSetSizes_Block( *A, 2, 2, 2, 2 );
#endif
    
    if(K) {MatBlockSetValue( *A, 0,0, K, SAME_NONZERO_PATTERN, INSERT_VALUES );}
    if(G) {MatBlockSetValue( *A, 0,1, G, SAME_NONZERO_PATTERN, INSERT_VALUES );}
    if(D) {MatBlockSetValue( *A, 1,0, D, SAME_NONZERO_PATTERN, INSERT_VALUES );}
    if(C) {MatBlockSetValue( *A, 1,1, C, SAME_NONZERO_PATTERN, INSERT_VALUES );}
    MatAssemblyBegin( *A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( *A, MAT_FINAL_ASSEMBLY );
    
    MatGetVecs( *A, x, b );
    
    if(u) {VecBlockSetValue( *x, 0, u, INSERT_VALUES );}
    if(p) {VecBlockSetValue( *x, 1, p, INSERT_VALUES );}
    VecAssemblyBegin( *x );
    VecAssemblyEnd( *x );
    
    if(f) {VecBlockSetValue( *b, 0, f, INSERT_VALUES );}
    if(h) {VecBlockSetValue( *b, 1, h, INSERT_VALUES );}
    VecAssemblyBegin( *b );
    VecAssemblyEnd( *b );
    
    
    PetscFunctionReturn(0);
}

void SBSNES_GetStokesOperators( 
		Stokes_SLE *stokesSLE,
		Mat *K,Mat *G,Mat *D,Mat *C,Mat *Smat,
		Vec *f,Vec *h,Vec *u,Vec *p )
{
	
    if(K) {if (stokesSLE->kStiffMat){      *K = SBSNES_GetPetscMatrix( stokesSLE->kStiffMat->matrix );     } else { *K = PETSC_NULL;}}
    if(G) {if (stokesSLE->gStiffMat){      *G = SBSNES_GetPetscMatrix( stokesSLE->gStiffMat->matrix );     } else { *G = PETSC_NULL;}}
    if(D) {if (stokesSLE->dStiffMat){      *D = SBSNES_GetPetscMatrix( stokesSLE->dStiffMat->matrix );     } else { *D = PETSC_NULL;}}
    if(C) {if (stokesSLE->cStiffMat){      *C = SBSNES_GetPetscMatrix( stokesSLE->cStiffMat->matrix );     } else { *C = PETSC_NULL;}}
	
	/* preconditioner */
	if(Smat){
	    *Smat = PETSC_NULL;
	    if( ((StokesBlockSNESInterface*)stokesSLE->solver)->preconditioner ) {
		StiffnessMatrix *preconditioner;
		
		preconditioner = ((StokesBlockSNESInterface*)stokesSLE->solver)->preconditioner;
		*Smat = SBSNES_GetPetscMatrix( preconditioner->matrix );
	    }
	}
	if(f){
	    *f = PETSC_NULL;
	    if (stokesSLE->fForceVec){      *f = SBSNES_GetPetscVector( stokesSLE->fForceVec->vector );     }
	}
	if(h){
	    *h = PETSC_NULL;
	    if (stokesSLE->hForceVec){      *h = SBSNES_GetPetscVector( stokesSLE->hForceVec->vector );     }
	}
	if(u){
	    *u = PETSC_NULL;
	    if (stokesSLE->uSolnVec){       *u = SBSNES_GetPetscVector( stokesSLE->uSolnVec->vector );      }
	}
	if(p){
	    *p = PETSC_NULL;
	    if (stokesSLE->pSolnVec){       *p = SBSNES_GetPetscVector( stokesSLE->pSolnVec->vector );      }
	}
}
void SBSNES_Assemble( StokesBlockSNESInterface* self, Bool nonLinear ) {
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
PetscErrorCode SBSNES_Check(Mat M, Vec V, const char *f, char* Mname, char* Vname){
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
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/* to-do: these two functions need to be somewhere common now like a 'Utils' directory */
void sb_writeMat(Mat K, char name[], char message[]){
    PetscViewer             mat_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
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
}
void sb_writeVec(Vec V, char name[], char message[]){
    PetscViewer             vec_view_file;
    char str[100];
    PetscTruth dump=PETSC_FALSE;
    PetscOptionsGetTruth( PETSC_NULL ,"-dump_matvec", &dump, 0 );
    if(dump){
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
}
PetscErrorCode SBSNES_FormResidual( SNES snes, Vec w, Vec F, void* _self )
{
   Mat K, G, D, C;
   Vec u,p,f,h, *subVecs;
   //Vec save_u, save_p;/* temporary vectors to hold current solution */
   double wallTime;// uNorm, pNorm, bNorm;
   static double totTime=0.0;
   static int timesCalled=0;
   PetscTruth useTimer, found, flg;// uDiff=PETSC_FALSE;
   //PetscErrorCode          ierr;
   StokesBlockSNESInterface* self = (StokesBlockSNESInterface*) _self;
   Stokes_SLE* SLE = (Stokes_SLE*)(self->st_sle);

   useTimer = PETSC_FALSE;
   PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_use_timer", &useTimer, &found );
   if(useTimer){
       wallTime = MPI_Wtime();
   }

   /* w is current solution being iterated on */
   /* u and p are the SLE solution vectors */
   /* so here we are setting the SLE vectors to the current solution */
   /* update solution onto nodes  happens in the assemble function */
   SBSNES_GetStokesOperators( self->st_sle, &K, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL,  PETSC_NULL, PETSC_NULL, &u, &p );
   /* SBSNES_Check(0,u, __func__, 0, "u"); */
   /* SBSNES_Check(0,p, __func__, 0, "p"); */

   VecBlockGetSubVectors( w, &subVecs );
   /* SBSNES_Check(0,subVecs[0], __func__, 0, "w[0]"); */
   /* SBSNES_Check(0,subVecs[1], __func__, 0, "w[1]"); */
   VecCopy( subVecs[0], u );
   VecCopy( subVecs[1], p );
   VecBlockRestoreSubVectors( w );
   
   /* need to pass norms back down through assembly to rheology */
#if !( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5) )
   SNESGetFunctionNorm(snes, &(SLE->fnorm));
#else
   Vec vec_func;
   SNESGetFunction(snes,&vec_func,0,0);
   VecNorm(vec_func,NORM_2, &(SLE->fnorm) );
#endif
   SLE->fnorm=sqrt(SLE->fnorm);/* more or less guessing what current norm might be -- snes->norm is the previous residual norm */
   MatNorm(K,NORM_1,&(SLE->knorm));

   SBSNES_GetStokesOperators( self->st_sle, &K, &G, &D, &C, PETSC_NULL, &f, &h, PETSC_NULL, PETSC_NULL );
   flg = PETSC_FALSE;
   PetscOptionsGetTruth(PETSC_NULL,"-snes_dump_residual",&flg,0);
   if(flg){
       sb_writeMat( K, "Koo", "Dumping original K matrix in FormResidual Function");
       sb_writeMat( D, "Doo", "Dumping D matrix");
       sb_writeMat( G, "Goo", "Dumping G matrix");
       sb_writeVec( u, "uoo", "Dumping current guessed u Vector");
       sb_writeVec( p, "poo", "Dumping current guessed p Vector");
       sb_writeVec( f, "foo", "Dumping original f Vector");
       sb_writeVec( h, "hoo", "Dumping original h Vector");       
   }

   SBSNES_Assemble( self, False );

   SBSNES_GetStokesOperators( self->st_sle, &K, &G, &D, &C, PETSC_NULL, &f, &h, PETSC_NULL, PETSC_NULL );
   flg = PETSC_FALSE;
   PetscOptionsGetTruth(PETSC_NULL,"-snes_dump_residual",&flg,0);
   if(flg){
       sb_writeMat( K, "Ko", "Dumping original K matrix in FormResidual Function");
       sb_writeMat( D, "Do", "Dumping D matrix");
       sb_writeMat( G, "Go", "Dumping G matrix");
       sb_writeVec( u, "uo", "Dumping current guessed u Vector");
       sb_writeVec( p, "po", "Dumping current guessed p Vector");
       sb_writeVec( f, "fo", "Dumping original f Vector");
       sb_writeVec( h, "ho", "Dumping original h Vector");       
   }
   VecBlockGetSubVectors( F, &subVecs );
   VecScale( f, -1.0 );
   VecScale( h, -1.0 );
   MatMultAdd( K, u, f, f ); /* K*u + f -> f */
   MatMultAdd( G, p, f, f );
   if( D ) MatMultAdd( D, u, h, h );
   else MatMultTransposeAdd( G, u, h, h );
   if( C ) MatMultAdd( C, p, h, h );

   VecBlockGetSubVectors( F, &subVecs );
   VecCopy( f, subVecs[0] );
   VecCopy( h, subVecs[1] );
   VecBlockRestoreSubVectors( F );

   VecAssemblyBegin( F );
   VecAssemblyEnd( F );

   flg = PETSC_FALSE;
   PetscOptionsGetTruth(PETSC_NULL,"-snes_dump_residual",&flg,0);
   if(flg){
       sb_writeVec( F, "F", "Dumping Jacobian Residual (rhs) Vector");
       sb_writeVec( f, "fj", "Dumping f (jacobian rhs) Vector");
       sb_writeVec( h, "hj", "Dumping h (jacobian rhs) Vector");
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
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5) )
PetscErrorCode SBSNES_FormJacobian(
   SNES snes, Vec x, Mat J, Mat Jp, void* _self )
#else
PetscErrorCode SBSNES_FormJacobian(
   SNES snes, Vec x, Mat* J, Mat* Jp, MatStructure* flag, void* _self )
#endif
{
    Mat K, G, D, C, Gt;
    Vec u, p, *subVecs;
    //Vec save_u, save_p;/* temporary vectors to hold current solution */
    double wallTime;
    static double totTime=0.0;
    static int timesCalled=0;
    PetscTruth useTimer, found, flg;// uDiff=PETSC_FALSE;
    //PetscErrorCode          ierr;
    //PetscInt j,start,end;
    //PetscScalar val, amp, norm;
    StokesBlockSNESInterface* self = (StokesBlockSNESInterface*) _self;
    UnderworldContext*  context = (UnderworldContext*)self->ctx;

    //FiniteElementContext_CalcNewDt( context );

   useTimer = PETSC_FALSE;
   PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_use_timer", &useTimer, &found );
   if(useTimer){
       wallTime = MPI_Wtime();
   }

   SBSNES_GetStokesOperators( self->st_sle, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL,
			     &u, &p );
   /* x is current solution being iterated on */
   /* u and p are the SLE solution vectors */
   /* SBSNES_Check(0,u, __func__, 0, "u"); */
   /* SBSNES_Check(0,p, __func__, 0, "p"); */
   VecBlockGetSubVectors( x, &subVecs );
   /* SBSNES_Check(0,subVecs[0], __func__, 0, "x[0]"); */
   /* SBSNES_Check(0,subVecs[1], __func__, 0, "x[1]"); */
   VecCopy( subVecs[0], u );
   VecCopy( subVecs[1], p );
   VecBlockRestoreSubVectors( x );

   flg=PETSC_FALSE;
   PetscOptionsGetTruth( PETSC_NULL ,"-experimental_snes_use_timer", &flg, &found );
   if(flg){
       
   }
   /* this is getting twice just for the moment */
   SystemLinearEquations_UpdateSolutionOntoNodes( self->st_sle, self->ctx );
   FiniteElementContext_CalcNewDt( context );
   
   SBSNES_Assemble( self, True );
   SBSNES_GetStokesOperators( self->st_sle, &K, &G, &D, &C, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL );
   /* if(amp > 1e-50 && it ==0){ */
   /*     //VecCopy(save_u, u); */
   /*     Stg_VecDestroy(&save_u); */
   /* } */

   /* create a symbolic Gt if no D */
   if( !D ) {
       MatTranspose( G, MAT_INITIAL_MATRIX, &Gt);
       self->DIsSym = PETSC_TRUE;
   }
   else {
       Gt = D;
       self->DIsSym = PETSC_FALSE;
   }

   flg = PETSC_FALSE;
   PetscOptionsGetTruth(PETSC_NULL,"-snes_dump_suboperators",&flg,0);
   if(flg){
       sb_writeMat( K, "K", "Dumping K Matrix");
       sb_writeMat( G, "G", "Dumping G Matrix");
       sb_writeMat( Gt, "Gt", "Dumping Gt Matrix");
       sb_writeMat( C, "C", "Dumping C Matrix");
   }

   SBSNES_FormBlockOperator( K,G,Gt,C, &J );

   if( !D ) { Stg_MatDestroy(&Gt ); }/* relinquish control of Gt from this function entirely to J matrix */

   if(useTimer){
       wallTime = MPI_Wtime() - wallTime;
       PetscPrintf( PETSC_COMM_WORLD, "  J assemble time = %f \n", wallTime );
       totTime += wallTime;
       timesCalled++;
       PetscPrintf( PETSC_COMM_WORLD, "  J assemble time = %f : Accumulated Total. Called %d times \n", totTime, timesCalled );
   }

   return 0;
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/* Sets up Solver to be a custom ksp (KSP_BSSCR) solve by default: */
/* requires PetscExt */
void _StokesBlockSNESInterface_Solve( void* solver, void* _stokesSLE ) {
        Stokes_SLE*  stokesSLE  = (Stokes_SLE*)_stokesSLE;
	StokesBlockSNESInterface* self    = (StokesBlockSNESInterface*)solver;

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
	//Stream*   errorStream = Journal_Register( Error_Type, (Name)StokesBlockSNESInterface_Type  );
	//Mat stokes_P;
	Mat stokes_A=0;
	Vec stokes_x;
	//Vec stokes_b;

	SNES snes;
	KSP stokes_ksp;
	PC  stokes_pc;
        //char name[100];
	
	PetscTruth sym, flg;
	static int called=0;

	SBSNES_GetStokesOperators( stokesSLE, &K,&G,&D,&C, &Smat, &f,&h, &u,&p );
	if(called==0){
	    PetscScalar *array;
	    PetscInt rstart, rend, k, i;
	    VecGetOwnershipRange(u,&rstart,&rend);
	    VecGetArray(u,&array);
	    k = 0;
	    for (i=rstart; i<rend; i++){
		array[k] = (0.01*rand()/(RAND_MAX+1.0));
		k++;
	    }
	    VecRestoreArray(u,&array);
	}
	called++;

        /* create a symbolic Gt */
	if( !D ) {
        MatTranspose( G, MAT_INITIAL_MATRIX, &Gt);
	    sym = PETSC_TRUE;
	    self->DIsSym = sym;
	}
	else {
	    Gt = D;
	    sym = PETSC_FALSE;
	    self->DIsSym = sym;
	}
	sb_writeMat( D, "Dstokes", "Dumping D matrix");
	SBSNES_FormBlockSystem( K,G,Gt,C, u,p, f,h, &stokes_A, &stokes_x, &F );

	SNESCreate( PETSC_COMM_WORLD, &snes );
	SNESSetJacobian( snes, stokes_A, stokes_A, SBSNES_FormJacobian, (void*)self );
	SNESSetFunction( snes, F, SBSNES_FormResidual, (void*)self );

	SNESGetKSP( snes, &stokes_ksp );
	KSPSetType( stokes_ksp, "bsscr" );/* i.e. making this the default solver : calls KSPCreate_XXX */
	KSPGetPC( stokes_ksp, &stokes_pc );
	PCSetType( stokes_pc, PCNONE );
	KSPSetFromOptions( stokes_ksp );

	/*
	   Doing this so the KSP Solver has access to the StgFEM Multigrid struct (PETScMGSolver).
	   As well as any custom stuff on the Stokes_SLE struct 
	*/
	if( stokes_ksp->data ){/* then ksp->data has been created in a KSpSetUp_XXX function */
	    /* testing for our KSP types that need the data that is on Solver... */
	    /* for the moment then, this function not completely agnostic about our KSPs */
	    flg=PETSC_FALSE;
	    PetscOptionsHasName(PETSC_NULL,"-use_petsc_ksp",&flg);
	    if (!flg) {
		((KSP_COMMON*)(stokes_ksp->data))->st_sle         = self->st_sle;
		((KSP_COMMON*)(stokes_ksp->data))->mg             = self->mg;
		((KSP_COMMON*)(stokes_ksp->data))->DIsSym         = self->DIsSym;
		((KSP_COMMON*)(stokes_ksp->data))->preconditioner = self->preconditioner;
	    }
	}
	//KSPSetInitialGuessNonzero( stokes_ksp, PETSC_FALSE );
	SNESSetFromOptions( snes );
	PetscReal abstol, rtol, bnorm, fnorm, hnorm;
	PetscReal stol;
	PetscInt maxit;
	PetscInt maxf;
	VecNorm(f,NORM_2,&fnorm);
	VecNorm(h,NORM_2,&hnorm);
	bnorm=sqrt(fnorm*fnorm+hnorm*hnorm);
	SNESGetTolerances( snes, &abstol, &rtol, &stol, &maxit, &maxf );
	abstol=rtol*bnorm*0.1;
	SNESSetTolerances( snes, abstol, rtol, stol, maxit, maxf );
	SNESSetFromOptions( snes ); /* so can still override the abstol that is set here */
	SNESGetTolerances( snes, &abstol, &rtol, &stol, &maxit, &maxf );

	/* Vec *subVecs; */
	/* PetscReal norm, amp, val; */
	/* PetscInt start, end, j; */
	/* amp=0.0; */
	/* PetscOptionsGetReal(PETSC_NULL,"-rand_amp",&amp,PETSC_NULL); */ 
	/* if(amp > 1e-50){ */
	/*     //VecDuplicate(u,&save_u); */
	/*     //VecCopy(u,save_u); */
	/*     //SNESGetFunctionNorm( snes, &norm); */
	/*     VecNorm(u, NORM_2, &norm); */
	/*     VecGetOwnershipRange(u,&start,&end); */
	/*     for(j=start;j<end;j++){ */
	/* 	if(norm < 1e-2){ */
	/* 	    val=norm*amp*(2.0*rand()/(RAND_MAX+1.0)-1.0); */
	/* 	} */
	/* 	else{ */
	/* 	    val=amp*(2.0*rand()/(RAND_MAX+1.0)-1.0); */
	/* 	} */
	/* 	VecSetValue( u, j, val, ADD_VALUES); */
	/*     } */
	/*     VecAssemblyBegin(u); */
	/*     VecAssemblyEnd(u); */
	/* } */
	    sb_writeVec( u, "Up", "Dumping prev Velocity Solution Vector");
	    sb_writeVec( p, "Pp", "Dumping prev Pressure Solution Vector");

	SNESSolve( snes, PETSC_NULL, stokes_x );
	flg = PETSC_FALSE;
	//PetscOptionsGetBool(PETSC_NULL,"-snes_dump_solution",&flg,PETSC_NULL);
	PetscOptionsGetTruth(PETSC_NULL,"-snes_dump_solution",&flg,0);
	if(flg){
	    sb_writeVec( u, "U", "Dumping Velocity Solution Vector");
	    sb_writeVec( p, "P", "Dumping Pressure Solution Vector");
	}

	Stg_SNESDestroy(&snes );

	Stg_MatDestroy(&stokes_A );
	Stg_VecDestroy(&stokes_x );
	Stg_VecDestroy(&F );

	if(!D){ Stg_MatDestroy(&Gt); }
}
#endif
