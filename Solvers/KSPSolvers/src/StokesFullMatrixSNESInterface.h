#ifndef __StokesFullMatrixSNESInterface_h__
#define __StokesFullMatrixSNESInterface_h__


//#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==2) )
	/** Textual name of this class */
	extern const Type StokesFullMatrixSNESInterface_Type;

	#define __StokesFullMatrixSNESInterface \
		/* General info */ \
		__SLE_Solver \
		/* Virtual info */ \
		StiffnessMatrix* preconditioner; \
		/* StokesFullMatrixSNESInterface info */ \
		Stokes_SLE *       st_sle;  \
                PETScMGSolver *    mg;  \
                int DIsSym; \
                UnderworldContext* ctx;

	struct StokesFullMatrixSNESInterface { __StokesFullMatrixSNESInterface };
	
	/* --- Constructors / Destructor --- */

	/** Constructor */
	void* _StokesFullMatrixSNESInterface_DefaultNew( Name name );
		
	/** Creation implementation / Virtual constructor */
	
	#define STOKESFULLMATRIXSNESINTERFACE_DEFARGS \
                SLE_SOLVER_DEFARGS

	#define STOKESFULLMATRIXSNESINTERFACE_PASSARGS \
                SLE_SOLVER_PASSARGS

	StokesFullMatrixSNESInterface* _StokesFullMatrixSNESInterface_New(  STOKESFULLMATRIXSNESINTERFACE_DEFARGS  );

	/** Class member variable initialisation */
	void _StokesFullMatrixSNESInterface_Init(
		StokesFullMatrixSNESInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg );
		

	/** Stg_Component_Build() implementations: allocates the 2 MatrixSolvers and additional Vectors */
	void _StokesFullMatrixSNESInterface_Build( void* solver, void* stokesSLE );
	
	void _StokesFullMatrixSNESInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data );
	
	void _StokesFullMatrixSNESInterface_Initialise( void* solver, void* stokesSLE ) ;
	
        /* void _StokesFullMatrixSNESInterface_Destroy( void* solver, void* data ); */

	void _StokesFullMatrixSNESInterface_SolverSetup( void* stokesSle, void* stokesSLE );

	void _StokesFullMatrixSNESInterface_Solve( void* solver, void* stokesSLE );
        PetscErrorCode SFMSNES_Check(Mat M, Vec V, const char *f, char* Mname, char* Vname);
//#endif	
#endif


