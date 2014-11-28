#ifndef __PICellerator_Common_Ppc_SphericalDepth_h__
#define __PICellerator_Common_Ppc_SphericalDepth_h__

/** Textual name of this class */
extern const Type Ppc_SphericalDepth_Type;

/** Ppc_SphericalDepth class contents */
#define __Ppc_SphericalDepth \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		FeMesh          *mesh; \
		Swarm           *ms; /* for passive tracer swarm */ \
      double          refHeight;

struct Ppc_SphericalDepth
{
   __Ppc_SphericalDepth
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_SphericalDepth* _Ppc_SphericalDepth_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_SphericalDepth_Delete( void* _self );
void _Ppc_SphericalDepth_Print( void* _self, Stream* stream );
void* _Ppc_SphericalDepth_DefaultNew( Name name );
void _Ppc_SphericalDepth_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_SphericalDepth_Build( void* _self, void* data );
void _Ppc_SphericalDepth_Initialise( void* _self, void* data );
void _Ppc_SphericalDepth_Execute( void* _self, void* data );
void _Ppc_SphericalDepth_Destroy( void* _self, void* data );

void _Ppc_SphericalDepth_ReferenceHeight( Ppc_SphericalDepth* self, void* data);
/* Public functions */
int _Ppc_SphericalDepth_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

