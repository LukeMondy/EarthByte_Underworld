#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "types.h"
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsnes.h>
#include <petscis.h> 
#include <petscviewer.h>

#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include "petsc-private/kspimpl.h"   /*I "petscksp.h" I*/
#else
  #include "private/kspimpl.h"   /*I "petscksp.h" I*/
#endif

//#include "ksptypes.h"
#include "ksp-register.h"
#include "StokesBlockKSPInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Macro for checking number integrity - i.e. checks if number is infinite or "not a number" */
#define SBKSP_isGoodNumber( number ) ( (! isnan( number ) ) && ( ! isinf( number ) ) )
#define SBKSP_GetPetscMatrix( matrix ) ( (Mat)(matrix) )
#define SBKSP_GetPetscVector( vector ) ( (Vec)(vector) )
#define SBKSP_GetPetscKSP( solver ) ( (KSP)(solver)  )

const Type StokesBlockKSPInterface_Type = "StokesBlockKSPInterface";

void* _StokesBlockKSPInterface_DefaultNew( Name name ) {
    SizeT                                              _sizeOfSelf = sizeof(StokesBlockKSPInterface);
    Type                                                      type = StokesBlockKSPInterface_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLE_Solver_Delete;
    Stg_Class_PrintFunction*                                _print = _SLE_Solver_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLE_Solver_Copy;
    Stg_Component_BuildFunction*                            _build = _StokesBlockKSPInterface_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _StokesBlockKSPInterface_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLE_Solver_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLE_Solver_Destroy;
    SLE_Solver_GetResidualFunc*                       _getResidual = NULL;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _StokesBlockKSPInterface_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _StokesBlockKSPInterface_AssignFromXML;
    SLE_Solver_SolverSetupFunction*                   _solverSetup = _StokesBlockKSPInterface_SolverSetup;
    SLE_Solver_SolveFunction*                               _solve = _StokesBlockKSPInterface_Solve;

    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
    return (void*) _StokesBlockKSPInterface_New( STOKESBLOCKKSPINTERFACE_PASSARGS );
}


/* Creation implementation / Virtual constructor */
/* Set up function pointers */
StokesBlockKSPInterface* _StokesBlockKSPInterface_New( STOKESBLOCKKSPINTERFACE_DEFARGS )
{    
    StokesBlockKSPInterface* self;
    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(StokesBlockKSPInterface) );

    self = (StokesBlockKSPInterface*) _SLE_Solver_New( SLE_SOLVER_PASSARGS );

    /* Virtual info */
    return self;
}

void _StokesBlockKSPInterface_Init( 
		StokesBlockKSPInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg,
        Name   filename, 
        char * string )
{
	self->preconditioner = preconditioner;
	self->st_sle = st_sle;
	self->mg     = mg;
    self->optionsFile = filename;
    self->optionsString = string;
}

void _StokesBlockKSPInterface_Build( void* solver, void* sle ) {/* it is the sle here being passed in*/
	StokesBlockKSPInterface*	self  = (StokesBlockKSPInterface*)solver;
	
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

void _StokesBlockKSPInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data ) {
	StokesBlockKSPInterface* self         = (StokesBlockKSPInterface*) solver;
	//double                  tolerance;
	//Iteration_Index         maxUzawaIterations, minUzawaIterations;
	StiffnessMatrix*        preconditioner;
	//Bool                    useAbsoluteTolerance;
	//Bool                    monitor;
	Stokes_SLE *            st_sle;
	PETScMGSolver *         mg;
	Name                filename;
	char* 		        string;

	_SLE_Solver_AssignFromXML( self, cf, data );

	filename = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"OptionsFile", "" );
	string = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"OptionsString", "" );

	printf("****************************************************************\n");
	printf("1. %s - Adding file %s to the options\n",self->name,filename);
	printf("2. %s - Adding string %s to the options\n",self->name,string);
	printf("****************************************************************\n");

	PetscOptionsInsertFile(PETSC_COMM_WORLD, filename, PETSC_FALSE);
	PetscOptionsInsertString(string);


	preconditioner = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Preconditioner", StiffnessMatrix, False, data  );
	st_sle  = Stg_ComponentFactory_ConstructByName( cf, (Name)"stokesEqn", Stokes_SLE, True, data  ); 
	mg      = Stg_ComponentFactory_ConstructByName( cf, (Name)"mgSolver", PETScMGSolver, False, data);

	_StokesBlockKSPInterface_Init( self, preconditioner, st_sle, mg, filename, string );

}

void _StokesBlockKSPInterface_Initialise( void* solver, void* stokesSLE ) {
	StokesBlockKSPInterface* self = (StokesBlockKSPInterface*) solver;
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
	//KSPRegisterTEST("Solvers/KSPSolvers/src");/* not sure if the path matters much here. everything still worked even though I had it wrong */
	//KSPRegisterSBKSP("Solvers/KSPSolvers/src");
	KSPRegisterAllKSP("Solvers/KSPSolvers/src");
}

/* SolverSetup */

void _StokesBlockKSPInterface_SolverSetup( void* solver, void* stokesSLE ) {
	StokesBlockKSPInterface* self = (StokesBlockKSPInterface*) solver;
	//Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Stream_UnIndentBranch( StgFEM_Debug );
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

void SBKSP_GetStokesOperators( 
		Stokes_SLE *stokesSLE,
		Mat *K,Mat *G,Mat *D,Mat *C,Mat *approxS,
		Vec *f,Vec *h,Vec *u,Vec *p )
{
	
	*K = *G = *D = *C = PETSC_NULL;
	if (stokesSLE->kStiffMat){      *K = SBKSP_GetPetscMatrix( stokesSLE->kStiffMat->matrix );     }
	if (stokesSLE->gStiffMat){      *G = SBKSP_GetPetscMatrix( stokesSLE->gStiffMat->matrix );     }
	if (stokesSLE->dStiffMat){      *D = SBKSP_GetPetscMatrix( stokesSLE->dStiffMat->matrix );     }
	if (stokesSLE->cStiffMat){      *C = SBKSP_GetPetscMatrix( stokesSLE->cStiffMat->matrix );     }
	
	/* preconditioner */
	*approxS = PETSC_NULL;
	if( ((StokesBlockKSPInterface*)stokesSLE->solver)->preconditioner ) {
		StiffnessMatrix *preconditioner;

		preconditioner = ((StokesBlockKSPInterface*)stokesSLE->solver)->preconditioner;
		*approxS = SBKSP_GetPetscMatrix( preconditioner->matrix );
	}	
	
	*f = *h = PETSC_NULL;
	if (stokesSLE->fForceVec){      *f = SBKSP_GetPetscVector( stokesSLE->fForceVec->vector );     }
	if (stokesSLE->hForceVec){      *h = SBKSP_GetPetscVector( stokesSLE->hForceVec->vector );     }
	
	*u = *p = PETSC_NULL;
	if (stokesSLE->uSolnVec){       *u = SBKSP_GetPetscVector( stokesSLE->uSolnVec->vector );      }
	if (stokesSLE->pSolnVec){       *p = SBKSP_GetPetscVector( stokesSLE->pSolnVec->vector );      }
	
}
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/* Sets up Solver to be a custom ksp (KSP_BSSCR) solve by default: */
/* requires PetscExt */
PetscErrorCode _StokesBlockKSPInterface_Solve( void* solver, void* _stokesSLE ) {
  Stokes_SLE*  stokesSLE  = (Stokes_SLE*)_stokesSLE;
  StokesBlockKSPInterface* Solver    = (StokesBlockKSPInterface*)solver;

  /* Create shortcuts to stuff needed on sle */
  Mat       K;
  Mat       G;
  Mat       Gt;
  Mat       D;
  Mat       C;
  Mat       approxS;
  Vec       u;
  Vec       p;
  Vec       f;
  Vec       h;

  //Stream*   errorStream = Journal_Register( Error_Type, (Name)StokesBlockKSPInterface_Type  );

  Mat stokes_P;
  Mat stokes_A;
  Vec stokes_x;
  Vec stokes_b;
  Mat a[2][2];
  Vec x[2];
  Vec b[2];

  KSP stokes_ksp;
  PC  stokes_pc;
  IS  isr[2];
  IS  isc[2];
  PetscInt rowsizeK, colsizeK, rowsizeG, colsizeG, i, first,step, globalstartK, globalstartG, globalrowsizeK, globalcolsizeG;
  PetscInt globalendK, globalendG, globalcolsizeK;;
  //char name[100];

  PetscTruth sym,flg;
  PetscErrorCode ierr;

  const PetscInt *indices;

  SBKSP_GetStokesOperators( stokesSLE, &K,&G,&D,&C, &approxS, &f,&h, &u,&p );

  /* create a symbolic Gt */
  if( !D ) {
    //MatCreateSymTrans( PETSC_COMM_WORLD, G, &Gt );
    //PetscPrintf( PETSC_COMM_WORLD, "\t* BUILDING TRANSPOSE of G!! \n");
    //MatCreateTranspose( G, &Gt);
    ierr = MatTranspose( G, MAT_INITIAL_MATRIX, &Gt);CHKERRQ(ierr);
    sym = PETSC_TRUE;
    Solver->DIsSym = sym;
  }
  else {
    Gt = D;
    sym = PETSC_FALSE;
    Solver->DIsSym = sym;
  }
#if 0
  PetscMPIInt rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  ierr = MatGetLocalSize(K, &rowsizeK, &colsizeK);CHKERRQ(ierr);
  ierr = MatGetLocalSize(Gt, &rowsizeG, &colsizeG);CHKERRQ(ierr);
  ierr = MatGetSize(K, &globalrowsizeK, &globalcolsizeK);CHKERRQ(ierr);
  //ierr = MatGetSize(G, NULL, &globalcolsizeG);CHKERRQ(ierr);
  printf("\t [%d]:K  matrix local sizes: %dx%d\n",rank,rowsizeK,colsizeK);
  printf("\t [%d]:Gt matrix local sizes: %dx%d\n",rank,rowsizeG,colsizeG);
  MatGetOwnershipRange(K, &globalstartK, &globalendK);
  MatGetOwnershipRange(Gt, &globalstartG, &globalendG);
  printf("\t [%d]:global row indices: K, Gt: %d-%d %d-%d\n",rank,globalstartK,globalendK-1, globalstartG, globalendG-1);
  /*                                                        length, start_pos */
  rowsizeK=globalendK-globalstartK;
  rowsizeG=globalendG-globalstartG;
  printf("\t [%d]:Global row sizes for K G: %d %d\n",rank,rowsizeK,rowsizeG);
  ierr = ISCreateStride(PetscObjectComm((PetscObject) K), rowsizeK, globalstartK, 1, &isr[0]);CHKERRQ(ierr);
  ierr = ISCreateStride(PetscObjectComm((PetscObject) K), rowsizeG, globalrowsizeK+globalstartG, 1, &isr[1]);CHKERRQ(ierr);
  ierr = MatGetSize(Gt, NULL, &globalcolsizeG);CHKERRQ(ierr);
  //colsizeK=globalcolsizeK;
  //colsizeG=globalcolsizeG;
  MatGetOwnershipRange(G, &globalstartG, &globalendG);
  ierr = MatGetLocalSize(G, &rowsizeG, &colsizeG);CHKERRQ(ierr);
  ierr = ISCreateStride(PetscObjectComm((PetscObject) K), colsizeK, 0, 1, &isc[0]);CHKERRQ(ierr);
  ierr = ISCreateStride(PetscObjectComm((PetscObject) K), colsizeG, globalcolsizeK, 1, &isc[1]);CHKERRQ(ierr);
    

  printf("\t [%d]:rowK colG: %d %d\n",rank,rowsizeK,colsizeG );

  ierr = ISGetIndices(isr[0],&indices);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:Printing indices directly\n",rank);CHKERRQ(ierr);
  for (i=0; i<rowsizeK; i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:%D\n",rank,indices[i]);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(isr[0],&indices);CHKERRQ(ierr);

  ierr = ISGetIndices(isr[1],&indices);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:Printing indices directly\n",rank);CHKERRQ(ierr);
  for (i=0; i<colsizeG; i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:%D\n",rank,indices[i]);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(isr[1],&indices);CHKERRQ(ierr);

  ierr = ISStrideGetInfo(isr[0],&first,&step);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:isr[0]: %D %D\n",rank,first, step);CHKERRQ(ierr);
  printf("\t [%d]:--isr[0]: %d %d\n",rank,first,step);
  ierr = ISStrideGetInfo(isr[1],&first,&step);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t [%d]:isr[1]: %D %D\n",rank,first, step);CHKERRQ(ierr);
  printf("\t [%d]:--isr[1]: %d %d\n",rank,first,step);
  
  printf("[%d]: ISR[0]--------------\n",rank);
  ierr = ISView(isr[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  printf("[%d]: ISR[1]--------------\n",rank);
  ierr = ISView(isr[1],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  printf("[%d]: ISC[0]--------------\n",rank);
  ierr = ISView(isc[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  printf("[%d]: ISC[1]--------------\n",rank);
  ierr = ISView(isc[1],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  a[0][0]=K;  a[0][1]=G;
  a[1][0]=Gt; a[1][1]=C;
  ierr = MatCreateNest(PetscObjectComm((PetscObject) K), 2, NULL, 2, NULL, (Mat *)a, &stokes_A);CHKERRQ(ierr);
  ierr = MatAssemblyBegin( stokes_A, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
  ierr = MatAssemblyEnd( stokes_A, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);



  x[0]=u;
  x[1]=p;
  ierr = VecCreateNest(PetscObjectComm((PetscObject) u), 2, NULL, x, &stokes_x);CHKERRQ(ierr);
  ierr = VecAssemblyBegin( stokes_x );CHKERRQ(ierr);
  ierr = VecAssemblyEnd( stokes_x);CHKERRQ(ierr);
    
  b[0]=f;
  b[1]=h;
  ierr = VecCreateNest(PetscObjectComm((PetscObject) f), 2, NULL, b, &stokes_b);CHKERRQ(ierr);
  ierr = VecAssemblyBegin( stokes_b );CHKERRQ(ierr);
  ierr = VecAssemblyEnd( stokes_b);CHKERRQ(ierr);
    
  if( approxS ) {
    a[0][0]=K;    a[0][1]=G;
    a[1][0]=NULL; a[1][1]=approxS;
    ierr = MatCreateNest(PetscObjectComm((PetscObject) K), 2, NULL, 2, NULL, (Mat *)a, &stokes_P);CHKERRQ(ierr);
    ierr = MatAssemblyBegin( stokes_P, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
    ierr = MatAssemblyEnd( stokes_P, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
  }
  else {
    stokes_P = stokes_A;
  }
    
  /* probably should make a Destroy function for these two */
  /* Update options from file and/or string here so we can change things on the fly */
  PetscOptionsInsertFile(PETSC_COMM_WORLD, Solver->optionsFile, PETSC_FALSE);
  PetscOptionsInsertString(Solver->optionsString);

  ierr = KSPCreate( PETSC_COMM_WORLD, &stokes_ksp );CHKERRQ(ierr);
  Stg_KSPSetOperators( stokes_ksp, stokes_A, stokes_P, SAME_NONZERO_PATTERN );
  ierr = KSPSetType( stokes_ksp, "bsscr" );/* i.e. making this the default solver : calls KSPCreate_XXX */CHKERRQ(ierr);
    
  ierr = KSPGetPC( stokes_ksp, &stokes_pc );CHKERRQ(ierr);
  ierr = PCSetType( stokes_pc, PCNONE );CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero( stokes_ksp, PETSC_TRUE );CHKERRQ(ierr);
  ierr = KSPSetFromOptions( stokes_ksp );CHKERRQ(ierr);
    

/* #if( PETSC_VERSION_MAJOR < 3 ) */
/* 	PCBlock_SetBlockType( stokes_pc, PC_BLOCK_UPPER ); */
/* #elif( PETSC_VERSION >= 3 ) */
/* 	PCBlockSetBlockType( stokes_pc, PC_BLOCK_UPPER ); */
/* #endif */

  /*
    Doing this so the KSP Solver has access to the StgFEM Multigrid struct (PETScMGSolver).
    As well as any custom stuff on the Stokes_SLE struct 
  */
  if( stokes_ksp->data ){/* then ksp->data has been created in a KSpSetUp_XXX function */
    /* testing for our KSP types that need the data that is on Solver... */
    /* for the moment then, this function not completely agnostic about our KSPs */
    //if(!strcmp("bsscr",stokes_ksp->type_name)){/* if is bsscr then set up the data on the ksp */
    flg=PETSC_FALSE;
    PetscOptionsHasName(PETSC_NULL,"-use_petsc_ksp",&flg);
    if (!flg) {
      ((KSP_COMMON*)(stokes_ksp->data))->st_sle         = Solver->st_sle;
      ((KSP_COMMON*)(stokes_ksp->data))->mg             = Solver->mg;
      ((KSP_COMMON*)(stokes_ksp->data))->DIsSym         = Solver->DIsSym;
      ((KSP_COMMON*)(stokes_ksp->data))->preconditioner = Solver->preconditioner;
    }
  }

  ierr = KSPSolve( stokes_ksp, stokes_b, stokes_x );CHKERRQ(ierr);
    
  Stg_KSPDestroy(&stokes_ksp );
  if( ((StokesBlockKSPInterface*)stokesSLE->solver)->preconditioner ){ Stg_MatDestroy(&stokes_P ); }

  Stg_MatDestroy(&stokes_A );

  Stg_VecDestroy(&stokes_x);
  Stg_VecDestroy(&stokes_b);

  if(!D){ Stg_MatDestroy(&Gt); }

  PetscFunctionReturn(0);
}

