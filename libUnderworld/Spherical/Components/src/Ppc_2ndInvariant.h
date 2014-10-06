#ifndef __PICellerator_Common_Ppc_2ndInvariant_h__
#define __PICellerator_Common_Ppc_2ndInvariant_h__

/** Textual name of this class */
extern const Type Ppc_2ndInvariant_Type;

/** Ppc_2ndInvariant class contents */
#define __Ppc_2ndInvariant \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int tensorTag;

struct Ppc_2ndInvariant
{
   __Ppc_2ndInvariant
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_2ndInvariant* _Ppc_2ndInvariant_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_2ndInvariant_Delete( void* _self );
void _Ppc_2ndInvariant_Print( void* _self, Stream* stream );
void* _Ppc_2ndInvariant_DefaultNew( Name name );
void _Ppc_2ndInvariant_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_2ndInvariant_Build( void* _self, void* data );
void _Ppc_2ndInvariant_Initialise( void* _self, void* data );
void _Ppc_2ndInvariant_Execute( void* _self, void* data );
void _Ppc_2ndInvariant_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_2ndInvariant_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

