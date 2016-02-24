#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Components.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_Random.h"


/* Textual name of this class */
const Type Ppc_Random_Type = "Ppc_Random";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Random* _Ppc_Random_New(  PPC_OPERATION_DEFARGS  )
{
	Ppc_Random* self;
	
	assert( _sizeOfSelf >= sizeof(Ppc_Random) );
	nameAllocationType = NON_GLOBAL;
	self = (Ppc_Random*) _Ppc_New(  PPC_PASSARGS  );	
	self->_get = _get;
	return self;
}


void* _Ppc_Random_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Ppc_Random);
	Type                                                      type = Ppc_Random_Type;
	Stg_Class_DeleteFunction*                              _delete = _Ppc_Random_Delete;
	Stg_Class_PrintFunction*                                _print = _Ppc_Random_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Random_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Ppc_Random_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Ppc_Random_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Random_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Ppc_Random_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Ppc_Random_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Random_Get;

	return (void*)_Ppc_Random_New(  PPC_OPERATION_PASSARGS  );
}


void _Ppc_Random_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
  Ppc_Random* self = (Ppc_Random*)_self;
  
  /* Construct parent */
  _Ppc_AssignFromXML( self, cf, data );

  self->minrange = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"MinRange", "0" );
  self->maxrange = PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"MaxRange", "1" );
  
}


void _Ppc_Random_Build( void* _self, void* data ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Build parent */
  _Ppc_Build( self, data );
}

void _Ppc_Random_Initialise( void* _self, void* data ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Initialize parent */
  _Ppc_Initialise( self, data );
}

void _Ppc_Random_Delete( void* _self ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Delete parent */
  _Ppc_Delete( self );
}

void _Ppc_Random_Print( void* _self, Stream* stream ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Print parent */
  _Ppc_Print( self, stream );
}

void _Ppc_Random_Execute( void* _self, void* data ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Execute parent */
  _Ppc_Execute( self, data );
}

void _Ppc_Random_Destroy( void* _self, void* data ) {
  Ppc_Random* self = (Ppc_Random*)_self;

  /* Destroy parent */
  _Ppc_Destroy( self, data );
}


/* 
 * Public functions 
 *
 */

int _Ppc_Random_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
   Ppc_Random* self = (Ppc_Random*) _self;
   int randnum = rand();
   double minr, maxr;
   int err;
   double ratio;


   err = PpcManager_Get( self->manager, lElement_I, particle, self->minrange, &minr );
   err = PpcManager_Get( self->manager, lElement_I, particle, self->maxrange, &maxr );
   //srand( 0 );

   ratio = (double)randnum/(double)RAND_MAX;

   result[0] = minr + (ratio) * (maxr - minr);

   return 0;
}
