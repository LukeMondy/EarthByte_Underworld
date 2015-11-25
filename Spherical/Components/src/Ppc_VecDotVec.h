#ifndef __PICellerator_Common_Ppc_VecDotVec_h__
#define __PICellerator_Common_Ppc_VecDotVec_h__

/** Textual name of this class */
extern const Type Ppc_VecDotVec_Type;

/** Ppc_VecDotVec class contents */
#define __Ppc_VecDotVec \
   /* Parent info */ \
   __Ppc \
   /* General data */\
   int vec1;         \
   int vec2;         \
   Bool transformv1;

struct Ppc_VecDotVec {
   __Ppc_VecDotVec
};


#ifndef ZERO
#define ZERO 0
#endif

#define PPC_CONSTANT_DEFARGS \
                PPC_DEFARGS

#define PPC_CONSTANT_PASSARGS \
                PPC_PASSARGS

Ppc_VecDotVec* _Ppc_VecDotVec_New(  PPC_CONSTANT_DEFARGS  );

void _Ppc_VecDotVec_Delete( void* _self );
void _Ppc_VecDotVec_Print( void* _self, Stream* stream );
void* _Ppc_VecDotVec_DefaultNew( Name name );
void _Ppc_VecDotVec_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data );
void _Ppc_VecDotVec_Build( void* _self, void* data );
void _Ppc_VecDotVec_Initialise( void* _self, void* data );
void _Ppc_VecDotVec_Execute( void* _self, void* data );
void _Ppc_VecDotVec_Destroy( void* _self, void* data );

/* Public functions */
int _Ppc_VecDotVec_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result );

#endif

