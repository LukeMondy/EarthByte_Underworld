#ifndef __PICellerator_Common_Ppc_PartialMelt_Limited_h__
#define __PICellerator_Common_Ppc_PartialMelt_Limited_h__

/** Textual name of this class */
extern const Type Ppc_PartialMelt_Limited_Type;

/** Ppc_PartialMelt_Limited class contents */
#define __Ppc_PartialMelt_Limited						 \
/* Parent info */										 \
  __Ppc													 \
  /* General data */									 \
  int             liquidusTag;              \
  int             solidusTag;               \
  int             temperatureTag;				 \
  int             MeltLimitTag;				 \
 
struct Ppc_PartialMelt_Limited
{
   __Ppc_PartialMelt_Limited
};

#ifndef ZERO
#define ZERO 0
#endif

#define PPC_LINEARDENSITY_DEFARGS \
                PPC_DEFARGS

#define PPC_LINEARDENSITY_PASSARGS \
                PPC_PASSARGS

Ppc_PartialMelt_Limited* _Ppc_PartialMelt_Limited_New(  PPC_LINEARDENSITY_DEFARGS  );

void _Ppc_PartialMelt_Limited_Delete( void* _self );
void _Ppc_PartialMelt_Limited_Print( void* _self, Stream* stream );
void* _Ppc_PartialMelt_Limited_DefaultNew( Name name );
void _Ppc_PartialMelt_Limited_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_PartialMelt_Limited_Build( void* _self, void* data );
void _Ppc_PartialMelt_Limited_Initialise( void* _self, void* data );
void _Ppc_PartialMelt_Limited_Execute( void* _self, void* data );
void _Ppc_PartialMelt_Limited_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_PartialMelt_Limited_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

