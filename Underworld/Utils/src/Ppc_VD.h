#ifndef __PICellerator_Common_Ppc_VD_h__
#define __PICellerator_Common_Ppc_VD_h__

/** Textual name of this class */
extern const Type Ppc_VD_Type;

/** Ppc_VD class contents */
#define __Ppc_VD \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int strainRateTag; \
		int viscosityTag;

struct Ppc_VD
{
   __Ppc_VD
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_VD* _Ppc_VD_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_VD_Delete( void* _self );
void _Ppc_VD_Print( void* _self, Stream* stream );
void* _Ppc_VD_DefaultNew( Name name );
void _Ppc_VD_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_VD_Build( void* _self, void* data );
void _Ppc_VD_Initialise( void* _self, void* data );
void _Ppc_VD_Execute( void* _self, void* data );
void _Ppc_VD_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_VD_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

