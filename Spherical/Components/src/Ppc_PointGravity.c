#include <mpi.h>
#include <string.h>
#include <math.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Components.h"

/* Textual name of this class */
const Type Ppc_PointGravity_Type = "Ppc_PointGravity";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_PointGravity* _Ppc_PointGravity_New(  PPC_POINTGRAVITY_DEFARGS  )
{
	Ppc_PointGravity* self;
	
	assert( _sizeOfSelf >= sizeof(Ppc_PointGravity) );
	nameAllocationType = NON_GLOBAL;
	self = (Ppc_PointGravity*) _Ppc_New(  PPC_PASSARGS  );	
	self->_get = _get;
	return self;
}


void* _Ppc_PointGravity_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Ppc_PointGravity);
	Type                                                      type = Ppc_PointGravity_Type;
	Stg_Class_DeleteFunction*                              _delete = _Ppc_PointGravity_Delete;
	Stg_Class_PrintFunction*                                _print = _Ppc_PointGravity_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_PointGravity_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Ppc_PointGravity_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Ppc_PointGravity_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Ppc_PointGravity_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Ppc_PointGravity_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Ppc_PointGravity_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_PointGravity_Get;

	return (void*)_Ppc_PointGravity_New(  PPC_POINTGRAVITY_PASSARGS  );
}


void _Ppc_PointGravity_Init( void* _self, int alphaTag, double x, double y, double z )
{
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

#if 0
#endif
	/* Sanity check */
	Journal_Firewall( alphaTag >= 0, 
							self->error_stream, "\n\n\n"
							"Error in configuration of Ppc_PointGravity %s\n"
							"Alpha not set"
							"\n\n\n", self->name );

	self->alphaTag = alphaTag;
	self->point[0] = x;
	self->point[1] = y;
	self->point[2] = z;
}


void _Ppc_PointGravity_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
	Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Construct parent */
	_Ppc_AssignFromXML( self, cf, data );

   _Ppc_PointGravity_Init( 
      self, 
      PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"Alpha", "" ), 
      Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"xcoord", 0.0 ),
      Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"ycoord", 0.0 ),
      Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"zcoord", 0.0 ) );
}


void _Ppc_PointGravity_Build( void* _self, void* data ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Build parent */
	_Ppc_Build( self, data );

}

void _Ppc_PointGravity_Initialise( void* _self, void* data ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Initialize parent */
	_Ppc_Initialise( self, data );
}

void _Ppc_PointGravity_Delete( void* _self ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Delete parent */
	_Ppc_Delete( self );
}

void _Ppc_PointGravity_Print( void* _self, Stream* stream ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Print parent */
	_Ppc_Print( self, stream );
}

void _Ppc_PointGravity_Execute( void* _self, void* data ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Execute parent */
	_Ppc_Execute( self, data );
}

void _Ppc_PointGravity_Destroy( void* _self, void* data ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*)_self;

	/* Destroy parent */
	_Ppc_Destroy( self, data );
}


Ppc_PointGravity* Ppc_PointGravity_New( Name name, Stg_Component* _self ) {
  Ppc_PointGravity* self = (Ppc_PointGravity*) _self;
  return self;
}

/* 
 * Public functions 
 *
 */


int _Ppc_PointGravity_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
   Ppc_PointGravity* self = (Ppc_PointGravity*) _self;
   MaterialPointsSwarm *ms=NULL;
   double vector[3], gravity[3], global[3], mag, alpha;
   int err;

   // get global coord
   FeMesh_CoordLocalToGlobal( self->manager->integrationSwarm->mesh, lElement_I, particle->xi, global );

   // calc vector to self->point
   vector[0] = self->point[0]-global[0];
   vector[1] = self->point[1]-global[1];
   vector[2] = self->point[2]-global[2];

   if (self->manager->integrationSwarm->dim == 2) vector[2]=0;

   mag = StGermain_VectorMagnitude( vector, self->manager->integrationSwarm->dim );

   // make a unit vector
   vector[0] = vector[0] / mag;
   vector[1] = vector[1] / mag;
   vector[2] = vector[2] / mag;

   assert( fabs(StGermain_VectorMagnitude( vector, self->manager->integrationSwarm->dim )-1) < 1e-8 );

  /* get alpha */
   err = PpcManager_Get( self->manager, lElement_I, particle, self->alphaTag, &alpha );
   assert( !err );

   /* TODO: Gravity doesn't need to be a static vector type now, this implementation should be changed so
      Gravity could be defined as a general ppc */ 
   err = PpcManager_GetGravity( self->manager, lElement_I, particle, gravity );
   assert( !err );
 
   /* return unit vector */
   result[0] = gravity[1] * alpha * vector[0];
   result[1] = gravity[1] * alpha * vector[1];
   result[2] = gravity[1] * alpha * vector[2];

  return 0;
}
