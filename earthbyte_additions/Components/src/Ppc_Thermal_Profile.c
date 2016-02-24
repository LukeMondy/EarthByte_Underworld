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
#include "Ppc_Thermal_Profile.h"


/* Textual name of this class */
const Type Ppc_Thermal_Profile_Type = "Ppc_Thermal_Profile";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Thermal_Profile* _Ppc_Thermal_Profile_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_Thermal_Profile* self;

   assert( _sizeOfSelf >= sizeof(Ppc_Thermal_Profile) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_Thermal_Profile*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_Thermal_Profile_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_Thermal_Profile);
   Type                                                      type = Ppc_Thermal_Profile_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_Thermal_Profile_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_Thermal_Profile_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Thermal_Profile_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_Thermal_Profile_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_Thermal_Profile_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Thermal_Profile_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_Thermal_Profile_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_Thermal_Profile_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Thermal_Profile_Get;

   return (void*)_Ppc_Thermal_Profile_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_Thermal_Profile_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );

   self->start_coord = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StartCoord", 0.0 );
   self->end_coord   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"EndCoord", 100000.0 );
   self->min_temp    = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MinTemp", 0.0 );
   self->max_temp    = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MaxTemp", 1000.0 );
   self->A 			 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"LinearCoefficient", 0.0 );
   self->B			 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"ExponentialCoefficient1", 0.0 );
   self->C 			 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"ExponentialCoefficient2", 0.0 );
   self->axis 		 = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"Axis", 1 );
   
   

}


void _Ppc_Thermal_Profile_Build( void* _self, void* data )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Build parent */
   _Ppc_Build( self, data );

}

void _Ppc_Thermal_Profile_Initialise( void* _self, void* data )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );

}

void _Ppc_Thermal_Profile_Delete( void* _self )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_Thermal_Profile_Print( void* _self, Stream* stream )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_Thermal_Profile_Execute( void* _self, void* data )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_Thermal_Profile_Destroy( void* _self, void* data )
{
   Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_Thermal_Profile_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
	Ppc_Thermal_Profile* self = (Ppc_Thermal_Profile*) _self;
	Coord coord;
	double start_coord, end_coord, min_temp, max_temp, A, B, C, temp;
	unsigned int vertaxis;
	
	start_coord = self->start_coord;
	end_coord = self->end_coord;
	min_temp = self->min_temp;
	max_temp = self->max_temp;
	A = self->A;
	B = self->B;
	C = self->C; 
	vertaxis = self->axis;
	
    FeMesh_CoordLocalToGlobal( self->manager->mesh, lElement_I, particle->xi, coord );

	if( coord[vertaxis] > start_coord )
		temp = min_temp;
	else if( coord[vertaxis] <= start_coord && coord[vertaxis] > end_coord )
		temp = min_temp + A*(start_coord-coord[vertaxis]) + B*(1-exp(-C*(start_coord-coord[vertaxis])));
	else
		temp = max_temp;

    if( temp > max_temp )
        temp = max_temp;
    if( temp < min_temp )
        temp = min_temp;
        
    result[0] = temp;
    
	return 0;
}
