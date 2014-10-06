#ifndef __PICellerator_Common_Ppc_Compression_h__
#define __PICellerator_Common_Ppc_Compression_h__

/** Textual name of this class */
extern const Type Ppc_Compression_Type;

/** Ppc_Compression class contents */
#define __Ppc_Compression \
		/* Parent info */ \
		__Ppc \
		/* General data */ \
		int tensorTag; \
		int densityTag; \
		int velocityTag; \
		int gradRhoTag; \
		int invType;

struct Ppc_Compression
{
   __Ppc_Compression
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_Compression* _Ppc_Compression_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_Compression_Delete( void* _self );
void _Ppc_Compression_Print( void* _self, Stream* stream );
void* _Ppc_Compression_DefaultNew( Name name );
void _Ppc_Compression_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_Compression_Build( void* _self, void* data );
void _Ppc_Compression_Initialise( void* _self, void* data );
void _Ppc_Compression_Execute( void* _self, void* data );
void _Ppc_Compression_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_Compression_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

