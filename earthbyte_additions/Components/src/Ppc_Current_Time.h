#ifndef __PICellerator_Common_Ppc_Current_Time_h__
#define __PICellerator_Common_Ppc_Current_Time_h__

/** Textual name of this class */
extern const Type Ppc_Current_Time_Type;

/** Ppc_Current_Time class contents */
#define __Ppc_Current_Time				\
  /* Parent info */					\
  __Ppc									\
  /* General data */					\


	struct Ppc_Current_Time { __Ppc_Current_Time };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_OPERATION_DEFARGS \
                PPC_DEFARGS

	#define PPC_OPERATION_PASSARGS \
                PPC_PASSARGS

	Ppc_Current_Time* _Ppc_Current_Time_New(  PPC_OPERATION_DEFARGS  );
	
	void _Ppc_Current_Time_Delete( void* _self );
	void _Ppc_Current_Time_Print( void* _self, Stream* stream );
	void* _Ppc_Current_Time_DefaultNew( Name name );
   void _Ppc_Current_Time_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_Current_Time_Build( void* _self, void* data );
	void _Ppc_Current_Time_Initialise( void* _self, void* data );
	void _Ppc_Current_Time_Execute( void* _self, void* data );
	void _Ppc_Current_Time_Destroy( void* _self, void* data );


   /* Public functions */
	int _Ppc_Current_Time_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );
	
#endif

