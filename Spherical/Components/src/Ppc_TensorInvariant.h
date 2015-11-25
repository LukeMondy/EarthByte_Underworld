#ifndef __PICellerator_Common_Ppc_TensorInvariant_h__
#define __PICellerator_Common_Ppc_TensorInvariant_h__

/** Textual name of this class */
extern const Type Ppc_TensorInvariant_Type;

/** Ppc_TensorInvariant class contents */
#define __Ppc_TensorInvariant \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int tensorTag; \
		int densityTag; \
		int velocityTag; \
		int gradRhoTag; \
		int invType;

struct Ppc_TensorInvariant
{
   __Ppc_TensorInvariant
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_TensorInvariant* _Ppc_TensorInvariant_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_TensorInvariant_Delete( void* _self );
void _Ppc_TensorInvariant_Print( void* _self, Stream* stream );
void* _Ppc_TensorInvariant_DefaultNew( Name name );
void _Ppc_TensorInvariant_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_TensorInvariant_Build( void* _self, void* data );
void _Ppc_TensorInvariant_Initialise( void* _self, void* data );
void _Ppc_TensorInvariant_Execute( void* _self, void* data );
void _Ppc_TensorInvariant_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_TensorInvariant_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

