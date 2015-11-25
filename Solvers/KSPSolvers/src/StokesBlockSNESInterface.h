#ifndef __StokesBlockSNESInterface_h__
#define __StokesBlockSNESInterface_h__

	/** Textual name of this class */
	extern const Type StokesBlockSNESInterface_Type;

	#define __StokesBlockSNESInterface \
		/* General info */ \
		__SLE_Solver \
		/* Virtual info */ \
		StiffnessMatrix* preconditioner; \
		/* StokesBlockSNESInterface info */ \
		Stokes_SLE *       st_sle;  \
                PETScMGSolver *    mg;  \
                int DIsSym; \
                UnderworldContext* ctx;

	struct StokesBlockSNESInterface { __StokesBlockSNESInterface };
	
	/* --- Constructors / Destructor --- */

	/** Constructor */
	void* _StokesBlockSNESInterface_DefaultNew( Name name );
		
	/** Creation implementation / Virtual constructor */
	
	#define STOKESBLOCKSNESINTERFACE_DEFARGS \
                SLE_SOLVER_DEFARGS

	#define STOKESBLOCKSNESINTERFACE_PASSARGS \
                SLE_SOLVER_PASSARGS

	StokesBlockSNESInterface* _StokesBlockSNESInterface_New(  STOKESBLOCKSNESINTERFACE_DEFARGS  );

	/** Class member variable initialisation */
	void _StokesBlockSNESInterface_Init(
		StokesBlockSNESInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg );
		

	/** Stg_Component_Build() implementations: allocates the 2 MatrixSolvers and additional Vectors */
	void _StokesBlockSNESInterface_Build( void* solver, void* stokesSLE );
	
	void _StokesBlockSNESInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data );
	
	void _StokesBlockSNESInterface_Initialise( void* solver, void* stokesSLE ) ;
	
        /* void _StokesBlockSNESInterface_Destroy( void* solver, void* data ); */

	void _StokesBlockSNESInterface_SolverSetup( void* stokesSle, void* stokesSLE );

	void _StokesBlockSNESInterface_Solve( void* solver, void* stokesSLE );
        PetscErrorCode SBSNES_Check(Mat M, Vec V, const char *f, char* Mname, char* Vname);
#endif


