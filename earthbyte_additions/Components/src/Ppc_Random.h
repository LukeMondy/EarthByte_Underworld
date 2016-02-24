#ifndef __PICellerator_Common_Ppc_Random_h__
#define __PICellerator_Common_Ppc_Random_h__

/** Textual name of this class */
extern const Type Ppc_Random_Type;

/** Ppc_Random class contents */
#define __Ppc_Random			   \
  /* Parent info */				   \
  __Ppc							   \
  /* General data */			   \
  double minrange;                  \
  double maxrange;                  \


	struct Ppc_Random { __Ppc_Random };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_OPERATION_DEFARGS \
                PPC_DEFARGS

	#define PPC_OPERATION_PASSARGS \
                PPC_PASSARGS

	Ppc_Random* _Ppc_Random_New(  PPC_OPERATION_DEFARGS  );
	
	void _Ppc_Random_Delete( void* _self );
	void _Ppc_Random_Print( void* _self, Stream* stream );
	void* _Ppc_Random_DefaultNew( Name name );
    void _Ppc_Random_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_Random_Build( void* _self, void* data );
	void _Ppc_Random_Initialise( void* _self, void* data );
	void _Ppc_Random_Execute( void* _self, void* data );
	void _Ppc_Random_Destroy( void* _self, void* data );


   /* Public functions */
	int _Ppc_Random_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );
	
#endif

