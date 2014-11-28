
#ifndef __Spherical_SphericalPeriodicAdvector_h__
#define __Spherical_SphericalPeriodicAdvector_h__
	
	/* Textual name of this class */
	extern const Type SphericalPeriodicAdvector_Type;

	#define __SphericalPeriodicAdvector \
		__PeriodicBoundariesManager

	struct SphericalPeriodicAdvector { __SphericalPeriodicAdvector };

	void* _SphericalPeriodicAdvector_DefaultNew( Name name );

	SphericalPeriodicAdvector* SphericalPeriodicAdvector_New( 
		Name						name,
		PICelleratorContext*	context,
		Mesh*						mesh, 
		Swarm*					swarm,
		Dictionary*				dictionary );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SPHERICALPERIODICADVECTOR_DEFARGS \
                PERIODICBOUNDARIESMANAGER_DEFARGS

	#define SPHERICALPERIODICADVECTOR_PASSARGS \
                PERIODICBOUNDARIESMANAGER_PASSARGS

	SphericalPeriodicAdvector* _SphericalPeriodicAdvector_New(  PERIODICBOUNDARIESMANAGER_DEFARGS  );

	void _SphericalPeriodicAdvector_Init( SphericalPeriodicAdvector* self ); 

	void _SphericalPeriodicAdvector_AssignFromXML( void* periodicBCsManager, Stg_ComponentFactory* cf, void* data );
	
	void _SphericalPeriodicAdvector_Delete( void* context );

	void _SphericalPeriodicAdvector_Print( void* context, Stream* stream );

	void* _SphericalPeriodicAdvector_Copy( void* periodicBCsManager, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _SphericalPeriodicAdvector_Build( void* periodicBCsManager, void* data );

	void _SphericalPeriodicAdvector_Initialise( void* periodicBCsManager, void* data );

	void _SphericalPeriodicAdvector_Execute( void* periodicBCsManager, void* data );

	void _SphericalPeriodicAdvector_Destroy( void* periodicBCsManager, void* data );

	void SphericalPeriodicAdvector_AddPeriodicBoundary( void* periodicBCsManager, Axis axis );

	void SphericalPeriodicAdvector_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I );

#endif

