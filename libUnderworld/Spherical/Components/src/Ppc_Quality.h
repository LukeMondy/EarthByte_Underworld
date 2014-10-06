#ifndef __PICellerator_Common_Ppc_Quality_h__
#define __PICellerator_Common_Ppc_Quality_h__

/** Textual name of this class */
extern const Type Ppc_Quality_Type;

/** Ppc_Quality class contents */
#define __Ppc_Quality \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int jac_conTag;

struct Ppc_Quality
{
   __Ppc_Quality
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_Quality* _Ppc_Quality_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_Quality_Delete( void* _self );
void _Ppc_Quality_Print( void* _self, Stream* stream );
void* _Ppc_Quality_DefaultNew( Name name );
void _Ppc_Quality_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_Quality_Build( void* _self, void* data );
void _Ppc_Quality_Initialise( void* _self, void* data );
void _Ppc_Quality_Execute( void* _self, void* data );
void _Ppc_Quality_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_Quality_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

