/**
  This Ppc evaluates:

   c * exp( x )

  where c & x are ppcs

**/

#ifndef __PICellerator_Common_Ppc_Exponential_h__
#define __PICellerator_Common_Ppc_Exponential_h__

/** Textual name of this class */
extern const Type Ppc_Exponential_Type;


/** Ppc_Exponential class contents */
#define __Ppc_Exponential \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int x_tag;         \
		int c_tag;         \

struct Ppc_Exponential
{
   __Ppc_Exponential
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_EXPONENTIAL_DEFARGS \
                PPC_DEFARGS

#define PPC_EXPONENTIAL_PASSARGS \
                PPC_PASSARGS

Ppc_Exponential* _Ppc_Exponential_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_Exponential_Delete( void* _self );
void _Ppc_Exponential_Print( void* _self, Stream* stream );
void* _Ppc_Exponential_DefaultNew( Name name );
void _Ppc_Exponential_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_Exponential_Build( void* _self, void* data );
void _Ppc_Exponential_Initialise( void* _self, void* data );
void _Ppc_Exponential_Execute( void* _self, void* data );
void _Ppc_Exponential_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_Exponential_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

