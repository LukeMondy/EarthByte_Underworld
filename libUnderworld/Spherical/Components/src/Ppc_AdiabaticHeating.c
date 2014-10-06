#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include <string.h>
#include <math.h>

#include "types.h"
#include "Ppc_AdiabaticHeating.h"


/* Textual name of this class */
const Type Ppc_AdiabaticHeating_Type = "Ppc_AdiabaticHeating";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_AdiabaticHeating* _Ppc_AdiabaticHeating_New(  PPC_LINEARDENSITY_DEFARGS  )
{
	Ppc_AdiabaticHeating* self;
	
	assert( _sizeOfSelf >= sizeof(Ppc_AdiabaticHeating) );
	nameAllocationType = NON_GLOBAL;
	self = (Ppc_AdiabaticHeating*) _Ppc_New(  PPC_PASSARGS  );	
	self->_get = _get;
	return self;
}


void* _Ppc_AdiabaticHeating_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Ppc_AdiabaticHeating);
	Type                                                      type = Ppc_AdiabaticHeating_Type;
	Stg_Class_DeleteFunction*                              _delete = _Ppc_AdiabaticHeating_Delete;
	Stg_Class_PrintFunction*                                _print = _Ppc_AdiabaticHeating_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_AdiabaticHeating_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Ppc_AdiabaticHeating_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Ppc_AdiabaticHeating_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Ppc_AdiabaticHeating_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Ppc_AdiabaticHeating_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Ppc_AdiabaticHeating_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_AdiabaticHeating_Get;

	return (void*)_Ppc_AdiabaticHeating_New(  PPC_LINEARDENSITY_PASSARGS  );
}


void _Ppc_AdiabaticHeating_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
	Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Construct parent */
	_Ppc_AssignFromXML( self, cf, data );

   self->coeffTag = PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "coeff", "" );
   self->gravityTag = PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Gravity", "" );
   self->temperatureTag = PpcManager_GetField( self->manager, cf, (Stg_Component*)self, "Temperature", True );
   self->velocityTag = PpcManager_GetField( self->manager, cf, (Stg_Component*)self, "Velocity", True );
   self->refTemperatureTag = PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "ReferenceTemperature", "" );
}


void _Ppc_AdiabaticHeating_Build( void* _self, void* data ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Build parent */
	_Ppc_Build( self, data );
}

void _Ppc_AdiabaticHeating_Initialise( void* _self, void* data ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Initialize parent */
	_Ppc_Initialise( self, data );
}

void _Ppc_AdiabaticHeating_Delete( void* _self ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Delete parent */
	_Ppc_Delete( self );
}

void _Ppc_AdiabaticHeating_Print( void* _self, Stream* stream ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Print parent */
	_Ppc_Print( self, stream );
}

void _Ppc_AdiabaticHeating_Execute( void* _self, void* data ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Execute parent */
	_Ppc_Execute( self, data );
}

void _Ppc_AdiabaticHeating_Destroy( void* _self, void* data ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*)_self;

	/* Destroy parent */
	_Ppc_Destroy( self, data );
}

/* 
 * Public functions 
 *
 */

int _Ppc_AdiabaticHeating_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
   Ppc_AdiabaticHeating* self = (Ppc_AdiabaticHeating*) _self;
   int err;

   double coeff, T, T_s;
   double grav[3], vel[3], v_dot_g;


   err = PpcManager_Get( self->manager, lElement_I, particle, self->velocityTag, vel );
   assert(!err);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->temperatureTag, &T );
   assert(!err);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->gravityTag, grav );
   assert(!err);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->coeffTag, &coeff );
   assert(!err);

   err = PpcManager_Get( self->manager, lElement_I, particle, self->refTemperatureTag, &T_s );
   assert(!err);

   /* take the unit vector of gravity and take its dot product with velocity
      Assume the vectors are both in cartesian definition */
   StGermain_VectorNormalise( grav, self->manager->integrationSwarm->dim );
   v_dot_g = StGermain_VectorDotProduct( vel, grav, self->manager->integrationSwarm->dim );

   *result = coeff * v_dot_g * (T+T_s);

   return 0;
}
