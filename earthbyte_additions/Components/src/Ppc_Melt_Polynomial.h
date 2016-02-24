#ifndef __PICellerator_Common_Ppc_Melt_Polynomial_h__
#define __PICellerator_Common_Ppc_Melt_Polynomial_h__

/** Textual name of this class */
extern const Type Ppc_Melt_Polynomial_Type;

/** Ppc_Melt_Polynomial class contents */
#define __Ppc_Melt_Polynomial						 \
/* Parent info */										 \
  __Ppc													 \
  /* General data */									 \
  int             t0Tag;              \
  int             t1Tag;               \
  int             t2Tag;				 \
  int             t3Tag;				 \
  int			  pressureTag; 			\
 
struct Ppc_Melt_Polynomial
{
   __Ppc_Melt_Polynomial
};

#ifndef ZERO
#define ZERO 0
#endif

#define PPC_LINEARDENSITY_DEFARGS \
                PPC_DEFARGS

#define PPC_LINEARDENSITY_PASSARGS \
                PPC_PASSARGS

Ppc_Melt_Polynomial* _Ppc_Melt_Polynomial_New(  PPC_LINEARDENSITY_DEFARGS  );

void _Ppc_Melt_Polynomial_Delete( void* _self );
void _Ppc_Melt_Polynomial_Print( void* _self, Stream* stream );
void* _Ppc_Melt_Polynomial_DefaultNew( Name name );
void _Ppc_Melt_Polynomial_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_Melt_Polynomial_Build( void* _self, void* data );
void _Ppc_Melt_Polynomial_Initialise( void* _self, void* data );
void _Ppc_Melt_Polynomial_Execute( void* _self, void* data );
void _Ppc_Melt_Polynomial_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_Melt_Polynomial_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

