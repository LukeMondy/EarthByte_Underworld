#ifndef __PICellerator_Common_Ppc_ElementLine_h__
#define __PICellerator_Common_Ppc_ElementLine_h__

/** Textual name of this class */
extern const Type Ppc_ElementLine_Type;

/** Ppc_ElementLine class contents */
#define __Ppc_ElementLine \
		/* Parent info */ \
		__Ppc \
		/* General data */ \

struct Ppc_ElementLine
{
   __Ppc_ElementLine
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_ElementLine* _Ppc_ElementLine_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_ElementLine_Delete( void* _self );
void _Ppc_ElementLine_Print( void* _self, Stream* stream );
void* _Ppc_ElementLine_DefaultNew( Name name );
void _Ppc_ElementLine_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_ElementLine_Build( void* _self, void* data );
void _Ppc_ElementLine_Initialise( void* _self, void* data );
void _Ppc_ElementLine_Execute( void* _self, void* data );
void _Ppc_ElementLine_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_ElementLine_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

