#ifndef __PICellerator_Common_Ppc_AdiabaticHeating_h__
#define __PICellerator_Common_Ppc_AdiabaticHeating_h__

/** Textual name of this class */
extern const Type Ppc_AdiabaticHeating_Type;

/** Ppc_AdiabaticHeating class contents */
#define __Ppc_AdiabaticHeating						 \
/* Parent info */										 \
  __Ppc													 \
  /* General data */									 \
  int             refTemperatureTag;		 \
  int             velocityTag;						 \
  int             gravityTag;						 \
  int             temperatureTag;				 \
  int             coeffTag;					 \

	struct Ppc_AdiabaticHeating { __Ppc_AdiabaticHeating };
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_LINEARDENSITY_DEFARGS \
                PPC_DEFARGS

	#define PPC_LINEARDENSITY_PASSARGS \
                PPC_PASSARGS

	Ppc_AdiabaticHeating* _Ppc_AdiabaticHeating_New(  PPC_LINEARDENSITY_DEFARGS  );
	
	void _Ppc_AdiabaticHeating_Delete( void* _self );
	void _Ppc_AdiabaticHeating_Print( void* _self, Stream* stream );
	void* _Ppc_AdiabaticHeating_DefaultNew( Name name );
   void _Ppc_AdiabaticHeating_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_AdiabaticHeating_Build( void* _self, void* data );
	void _Ppc_AdiabaticHeating_Initialise( void* _self, void* data );
	void _Ppc_AdiabaticHeating_Execute( void* _self, void* data );
	void _Ppc_AdiabaticHeating_Destroy( void* _self, void* data );

   /* Public functions */
   int _Ppc_AdiabaticHeating_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

