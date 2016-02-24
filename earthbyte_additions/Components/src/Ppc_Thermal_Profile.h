#ifndef __PICellerator_Common_Ppc_Thermal_Profile_h__
#define __PICellerator_Common_Ppc_Thermal_Profile_h__

/** Textual name of this class */
extern const Type Ppc_Thermal_Profile_Type;

/** Ppc_Thermal_Profile class contents */
#define __Ppc_Thermal_Profile				\
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		double start_coord; \
		double end_coord; \
		double min_temp; \
		double max_temp; \
		double A; \
		double B; \
		double C; \
        unsigned axis; 

	struct Ppc_Thermal_Profile { __Ppc_Thermal_Profile };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_OPERATION_DEFARGS \
                PPC_DEFARGS

	#define PPC_OPERATION_PASSARGS \
                PPC_PASSARGS

	Ppc_Thermal_Profile* _Ppc_Thermal_Profile_New(  PPC_OPERATION_DEFARGS  );
	
	void _Ppc_Thermal_Profile_Delete( void* _self );
	void _Ppc_Thermal_Profile_Print( void* _self, Stream* stream );
	void* _Ppc_Thermal_Profile_DefaultNew( Name name );
   void _Ppc_Thermal_Profile_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_Thermal_Profile_Build( void* _self, void* data );
	void _Ppc_Thermal_Profile_Initialise( void* _self, void* data );
	void _Ppc_Thermal_Profile_Execute( void* _self, void* data );
	void _Ppc_Thermal_Profile_Destroy( void* _self, void* data );


   /* Public functions */
	int _Ppc_Thermal_Profile_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );
	
#endif

