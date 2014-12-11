#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_VecDotVec.h"


/* Textual name of this class */
const Type Ppc_VecDotVec_Type = "Ppc_VecDotVec";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_VecDotVec* _Ppc_VecDotVec_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_VecDotVec* self;

   assert( _sizeOfSelf >= sizeof(Ppc_VecDotVec) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_VecDotVec*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_VecDotVec_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_VecDotVec);
   Type                                                      type = Ppc_VecDotVec_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_VecDotVec_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_VecDotVec_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_VecDotVec_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_VecDotVec_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_VecDotVec_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_VecDotVec_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_VecDotVec_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_VecDotVec_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_VecDotVec_Get;

   return (void*)_Ppc_VecDotVec_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_VecDotVec_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );
   self->vec1 = PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Vector1", "" );
   self->vec2 = PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Vector2", "" );

   self->transformv1 = Stg_ComponentFactory_GetBool( cf, (Name)self->name, "TransformVec1", False );
}


void _Ppc_VecDotVec_Build( void* _self, void* data )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Build parent */
   _Ppc_Build( self, data );
}

void _Ppc_VecDotVec_Initialise( void* _self, void* data )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_VecDotVec_Delete( void* _self )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_VecDotVec_Print( void* _self, Stream* stream )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_VecDotVec_Execute( void* _self, void* data )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_VecDotVec_Destroy( void* _self, void* data )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_VecDotVec_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_VecDotVec* self = (Ppc_VecDotVec*) _self;
   double dot, vec1[3], vec2[3];
   int err;

   memset(vec1, 0, sizeof(double)*3);
   memset(vec2, 0, sizeof(double)*3);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->vec1, &(vec1[0]) );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->vec2, &(vec2[0]) );

   if( self->transformv1 ) {
      double xyz[3], rtp[3];
      memcpy( rtp, vec1, sizeof(double)*3 );

      FeMesh_CoordLocalToGlobal( self->manager->mesh, lElement_I, particle->xi, xyz );
      Spherical_VectorRTP2XYZ( rtp, xyz, self->manager->integrationSwarm->dim, vec1 ) ;
   }

   dot = vec1[0]*vec2[0] +
         vec1[1]*vec2[1] +
         vec1[2]*vec2[2] ;

   *result = dot;

   return 0;
}
