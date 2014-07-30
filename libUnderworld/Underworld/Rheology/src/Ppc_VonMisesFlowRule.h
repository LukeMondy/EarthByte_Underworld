#ifndef __PICellerator_Common_Ppc_VonMisesFlowRule_h__
#define __PICellerator_Common_Ppc_VonMisesFlowRule_h__

/** Textual name of this class */
extern const Type Ppc_VonMisesFlowRule_Type;

/** Ppc_VonMisesFlowRule class contents */
#define __Ppc_VonMisesFlowRule				\
  /* Parent info */					\
  __Ppc									\
  /* General data */					\
  SystemLinearEquations *sle; /* to know nonlinearity */\
  Index         hasYieldedIndex; /* index on material swarm extension */\
  SwarmVariable* hasYieldedVariable; /* pointer to swarm variable */\
  int           viscosityTag; /* latent heat of fusion */\
  int           strainRateInvTag; \
  int           yieldCritTag; \


	struct Ppc_VonMisesFlowRule { __Ppc_VonMisesFlowRule };

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PPC_OPERATION_DEFARGS \
                PPC_DEFARGS

	#define PPC_OPERATION_PASSARGS \
                PPC_PASSARGS

	Ppc_VonMisesFlowRule* _Ppc_VonMisesFlowRule_New(  PPC_OPERATION_DEFARGS  );
	
	void _Ppc_VonMisesFlowRule_Delete( void* _self );
	void _Ppc_VonMisesFlowRule_Print( void* _self, Stream* stream );
	void* _Ppc_VonMisesFlowRule_DefaultNew( Name name );
   void _Ppc_VonMisesFlowRule_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_VonMisesFlowRule_Build( void* _self, void* data );
	void _Ppc_VonMisesFlowRule_Initialise( void* _self, void* data );
	void _Ppc_VonMisesFlowRule_Execute( void* _self, void* data );
	void _Ppc_VonMisesFlowRule_Destroy( void* _self, void* data );


   /* Public functions */
	int _Ppc_VonMisesFlowRule_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );
	
#endif

