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
#include "Ppc_TensorInvariant.h"


/* Textual name of this class */
const Type Ppc_TensorInvariant_Type = "Ppc_TensorInvariant";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_TensorInvariant* _Ppc_TensorInvariant_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_TensorInvariant* self;

   assert( _sizeOfSelf >= sizeof(Ppc_TensorInvariant) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_TensorInvariant*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_TensorInvariant_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_TensorInvariant);
   Type                                                      type = Ppc_TensorInvariant_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_TensorInvariant_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_TensorInvariant_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_TensorInvariant_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_TensorInvariant_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_TensorInvariant_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_TensorInvariant_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_TensorInvariant_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_TensorInvariant_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_TensorInvariant_Get;

   return (void*)_Ppc_TensorInvariant_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_TensorInvariant_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;
   Bool found;

   /* The PpcManager */
	self->manager = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
	if( !self->manager  )
		self->manager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data );
   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );
   self->tensorTag = PpcManager_GetField( self->manager, cf, (Stg_Component*)self, "Tensor", True );

   found = Stg_ComponentFactory_TryInt( cf, self->name, (Dictionary_Entry_Key)"Invariant", &self->invType );
   Journal_Firewall( found, global_error_stream, "Error in %s:\nCan't find\n <param name=\"Invariant\"> </param>\n in \'%s\' definition. Add the order of invariant number you want\n", __func__, self->name); 
   if( self->invType > 2 || self->invType < 1 ) {
      Journal_Firewall( 0, global_error_stream,
         "Error in \'%s\' input: Only \'Invariant\' 1 or 2 can be calculated. You have %d\n", self->name, self->invType);
   }
}


void _Ppc_TensorInvariant_Build( void* _self, void* data )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Build parent */
   _Ppc_Build( self, data );
}

void _Ppc_TensorInvariant_Initialise( void* _self, void* data )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );
}

void _Ppc_TensorInvariant_Delete( void* _self )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_TensorInvariant_Print( void* _self, Stream* stream )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_TensorInvariant_Execute( void* _self, void* data )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_TensorInvariant_Destroy( void* _self, void* data )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

int _Ppc_TensorInvariant_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_TensorInvariant* self = (Ppc_TensorInvariant*) _self;
   int err;
   double inv, tensor[6];

   err = PpcManager_Get( self->manager, lElement_I, particle, self->tensorTag, tensor );
   assert(!err);

   switch( self->invType ) {
      case 1:
         SymmetricTensor_GetTrace( tensor, Mesh_GetDimSize(self->manager->mesh), &inv );
         break;
      case 2:
         inv = SymmetricTensor_2ndInvariant( tensor, Mesh_GetDimSize(self->manager->mesh) );
         break;
   }
   *result = inv;
   return 0;
}
