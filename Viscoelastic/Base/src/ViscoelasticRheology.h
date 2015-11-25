

#ifndef __Experimental_ViscoelasticRheology_h__
#define __Experimental_ViscoelasticRheology_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type ViscoelasticRheology_Type;
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __ViscoelasticRheology \
		/* Parent info */ \
		__Rheology \
		/* Virtual functions go here */ \
		/* Other Info */\
		double					elasticTimeStep;    \
		double					mu;					\
		ViscoelasticForceTerm* 	veForceTerm;					
	

	struct ViscoelasticRheology { __ViscoelasticRheology };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define VISCOELASTICRHEOLOGY_DEFARGS \
                RHEOLOGY_DEFARGS

	#define VISCOELASTICRHEOLOGY_PASSARGS \
                RHEOLOGY_PASSARGS

	ViscoelasticRheology* _ViscoelasticRheology_New(  VISCOELASTICRHEOLOGY_DEFARGS  );
	
	/* 'Stg_Component' implementations */
	void* _ViscoelasticRheology_DefaultNew( Name name ) ;
	void _ViscoelasticRheology_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
	void _ViscoelasticRheology_Initialise( void* rheology, void* data );

	void _ViscoelasticRheology_ModifyConstitutiveMatrix( 
			void*                                              rheology, 
			ConstitutiveMatrix*                                constitutiveMatrix,
			MaterialPointsSwarm*                               materialPointsSwarm,
			Element_LocalIndex                                 lElement_I,
			MaterialPoint*                                     materialPoint,
			Coord                                              xi ) ;

	double _ViscoelasticRheology_Timestep( ViscoelasticRheology* self, void* context );
	
#endif

