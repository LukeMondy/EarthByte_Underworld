#ifndef __PICellerator_Common_Ppc_Curie_Condition_h__
#define __PICellerator_Common_Ppc_Curie_Condition_h__

/** Textual name of this class */
extern const Type Ppc_Curie_Condition_Type;

/** Ppc_Curie_Condition class contents */
#define __Ppc_Curie_Condition				\
  /* Parent info */								\
  __Ppc											\
  /* General data */							\
  int               field;						\
  int               value;						\
  int 				magorient;      	        \
  int 				storage;			     	\
  int               initialOrientation;         \
  
	struct Ppc_Curie_Condition { __Ppc_Curie_Condition };

	Ppc_Curie_Condition* Ppc_Curie_Condition_New( Name name, Stg_Component* _self );
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define Ppc_Curie_Condition_DEFARGS \
                PPC_DEFARGS

	#define Ppc_Curie_Condition_PASSARGS \
                PPC_PASSARGS

	Ppc_Curie_Condition* _Ppc_Curie_Condition_New(  Ppc_Curie_Condition_DEFARGS  );
	
	void _Ppc_Curie_Condition_Delete( void* _self );
	void _Ppc_Curie_Condition_Print( void* _self, Stream* stream );
	void* _Ppc_Curie_Condition_DefaultNew( Name name );
   void _Ppc_Curie_Condition_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
	void _Ppc_Curie_Condition_Build( void* _self, void* data );
	void _Ppc_Curie_Condition_Initialise( void* _self, void* data );
	void _Ppc_Curie_Condition_Execute( void* _self, void* data );
	void _Ppc_Curie_Condition_Destroy( void* _self, void* data );

   /* Public functions */
   int _Ppc_Curie_Condition_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

