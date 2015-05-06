#ifndef __StokesBlockKSPInterface_h__
#define __StokesBlockKSPInterface_h__

	/** Textual name of this class */
	extern const Type StokesBlockKSPInterface_Type;

        /** if too much stuff ends up on this struct then possibly just
        extend the BSSCR struct instead at some point */
        #define __StokesBlockKSPInterface  \
                /* General info */ \
                __SLE_Solver	   \
		/* Virtual info */ \
		StiffnessMatrix* preconditioner; \
		/* StokesBlockKSPInterface info */ \
		Stokes_SLE    *    st_sle;  \
        PETScMGSolver *    mg;      \
        int DIsSym;                                              \
        /* approxS often same as preconditioner above */         \
		Mat K2, M, approxS, S;                                   \
        /* S1 = 1/sqrt(diag(K)); S2 = scaling for pressure */    \
		Vec S1,S2;                                               \
        /* file and string for petsc options */                  \
        Name optionsFile;                                        \
        char * optionsString;

	struct StokesBlockKSPInterface { __StokesBlockKSPInterface };

        #define __KSP_COMMON \
                Stokes_SLE     *   st_sle; \
                PETScMGSolver  *   mg;	   \
                PetscTruth         DIsSym; \
                StiffnessMatrix*   preconditioner;

        struct KSP_COMMON { __KSP_COMMON };
        typedef struct KSP_COMMON KSP_COMMON;


	/* --- Constructors / Destructor --- */

	/** Constructor */
	void* _StokesBlockKSPInterface_DefaultNew( Name name );
		
	/** Creation implementation / Virtual constructor */
	
	#define STOKESBLOCKKSPINTERFACE_DEFARGS \
                SLE_SOLVER_DEFARGS

	#define STOKESBLOCKKSPINTERFACE_PASSARGS \
                SLE_SOLVER_PASSARGS

	StokesBlockKSPInterface* _StokesBlockKSPInterface_New(  STOKESBLOCKKSPINTERFACE_DEFARGS  );

	/** Class member variable initialisation */
	void _StokesBlockKSPInterface_Init(
		StokesBlockKSPInterface*      self,
		StiffnessMatrix*   preconditioner,
		Stokes_SLE *       st_sle,
		PETScMGSolver *    mg,
        Name   filename, 
        char * string );
		

	/** Stg_Component_Build() implementations: allocates the 2 MatrixSolvers and additional Vectors */
	void _StokesBlockKSPInterface_Build( void* solver, void* stokesSLE );
	
	void _StokesBlockKSPInterface_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data );
	
	void _StokesBlockKSPInterface_Initialise( void* solver, void* stokesSLE ) ;
	
        /* void _StokesBlockKSPInterface_Destroy( void* solver, void* data ); */

	void _StokesBlockKSPInterface_SolverSetup( void* stokesSle, void* stokesSLE );

	PetscErrorCode _StokesBlockKSPInterface_Solve( void* solver, void* stokesSLE );

        void SBKSP_GetStokesOperators( 
		Stokes_SLE *stokesSLE,
		Mat *K,Mat *G,Mat *D,Mat *C,Mat *approxS,
		Vec *f,Vec *h,Vec *u,Vec *p );

#endif
