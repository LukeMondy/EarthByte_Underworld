#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Components/Components.h>

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_2ndInvariant.h"


/* Textual name of this class */
const Type Ppc_2ndInvariant_Type = "Ppc_2ndInvariant";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_2ndInvariant* _Ppc_2ndInvariant_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_2ndInvariant* self;

   assert( _sizeOfSelf >= sizeof(Ppc_2ndInvariant) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_2ndInvariant*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_2ndInvariant_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_2ndInvariant);
   Type                                                      type = Ppc_2ndInvariant_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_2ndInvariant_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_2ndInvariant_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_2ndInvariant_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_2ndInvariant_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_2ndInvariant_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_2ndInvariant_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_2ndInvariant_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_2ndInvariant_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_2ndInvariant_Get;

   return (void*)_Ppc_2ndInvariant_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_2ndInvariant_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* The PpcManager */
	self->manager = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
	if( !self->manager  )
		self->manager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data );
   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );
	self->tensorTag = PpcManager_GetField( self->manager, cf, (Stg_Component*)self, "Tensor", True );
}


void _Ppc_2ndInvariant_Build( void* _self, void* data )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Build parent */
   _Ppc_Build( self, data );
}

void _Ppc_2ndInvariant_Initialise( void* _self, void* data )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_2ndInvariant_Delete( void* _self )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_2ndInvariant_Print( void* _self, Stream* stream )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_2ndInvariant_Execute( void* _self, void* data )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_2ndInvariant_Destroy( void* _self, void* data )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_2ndInvariant_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_2ndInvariant* self = (Ppc_2ndInvariant*) _self;
   int err;
   double inv, tensor[6];

   err = PpcManager_Get( self->manager, lElement_I, particle, self->tensorTag, tensor );
   assert(!err);

   /* 
TODO:
- bit dodgy about the manager and mesh bit
- need a way to check the rank of the 'tensorTag' data, maybe itsn't a 2nd order tensor
*/
   inv = SymmetricTensor_2ndInvariant( tensor, Mesh_GetDimSize(self->manager->mesh) );

   *result = inv;

   return 0;
}
