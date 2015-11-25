#ifndef __Spherical_Components_ForceVector_h__
#define __Spherical_Components_ForceVector_h__
	
	
	/* Textual name of this class */
	extern const Type SphericalForceVector_Type;
	
	/* StiffnessMatrix information */
	#define __SphericalForceVector  \
		/* General info */ \
		__ForceVector \
		\
		/* Virtual info */ \
	  int totalDofsThisElement; \
		double* rotMat; /* storage for rotation matrix */ \
		double* tmpMat; /* storage for rotation matrix calculation, DGEMM needs a third matrix*/ \
	
	struct SphericalForceVector { __SphericalForceVector };
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SPHERICALFORCEVECTOR_DEFARGS \
                FORCEVECTOR_DEFARGS

	#define SPHERICALFORCEVECTOR_PASSARGS \
                FORCEVECTOR_PASSARGS

	void* _SphericalForceVector_DefaultNew( Name name );

	SphericalForceVector* _SphericalForceVector_New(  SPHERICALFORCEVECTOR_DEFARGS  );

	void _SphericalForceVector_AssignFromXML( void* forceVector, Stg_ComponentFactory* cf, void* data );

	void _SphericalForceVector_Build( void* forceVector, void* data );

	void _SphericalForceVector_Destroy( void* forceVector, void* data );

	void SphericalForceVector_GlobalAssembly_General( void* forceVector ) ;

	void SphericalForceVector_AssembleElement( void* forceVector, Element_LocalIndex element_lI, double* elForceVecToAdd ) ;

#endif

