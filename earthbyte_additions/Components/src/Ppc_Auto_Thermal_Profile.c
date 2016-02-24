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
#include "Ppc_Auto_Thermal_Profile.h"


/* Textual name of this class */
const Type Ppc_Auto_Thermal_Profile_Type = "Ppc_Auto_Thermal_Profile";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_Auto_Thermal_Profile* _Ppc_Auto_Thermal_Profile_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_Auto_Thermal_Profile* self;

   assert( _sizeOfSelf >= sizeof(Ppc_Auto_Thermal_Profile) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_Auto_Thermal_Profile*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_Auto_Thermal_Profile_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_Auto_Thermal_Profile);
   Type                                                      type = Ppc_Auto_Thermal_Profile_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_Auto_Thermal_Profile_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_Auto_Thermal_Profile_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_Auto_Thermal_Profile_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_Auto_Thermal_Profile_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_Auto_Thermal_Profile_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_Auto_Thermal_Profile_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_Auto_Thermal_Profile_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_Auto_Thermal_Profile_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_Auto_Thermal_Profile_Get;

   return (void*)_Ppc_Auto_Thermal_Profile_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_Auto_Thermal_Profile_Init( void* _self,
        int     diffusivityTag,
        int     densityTag,
        int     cpTag,
        int     radiogenicTag,
        double  start_coord,
        double  end_coord,
        double  min_temp,
        double  max_temp,
        double  axis ) {
    Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;
    
    self->diffusivityTag = diffusivityTag;
    self->densityTag = densityTag;
    self->cpTag = cpTag;
    self->radiogenicTag = radiogenicTag;  
    
    self->start_coord = start_coord;
    self->end_coord = end_coord;
    self->min_temp = min_temp;
    self->max_temp = max_temp;
    self->axis = axis;
 
}


void _Ppc_Auto_Thermal_Profile_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );
   
   _Ppc_Auto_Thermal_Profile_Init(
        self,
        PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Diffusivity", ""),
        PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Density", "" ),
        PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "Cp", "" ),
        PpcManager_GetPpcFromDict( self->manager, cf, (Name)self->name, "RadiogenicHeat", "" ),
        Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StartCoord", 0.0 ),
        Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"EndCoord", 100000.0 ),
        Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MinTemp", 0.0 ),
        Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MaxTemp", 1000.0 ),
        Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"Axis", 1 ) );

}


void _Ppc_Auto_Thermal_Profile_Build( void* _self, void* data )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   Stg_Component_Build( self->manager, data, False );
   
   /* Build parent */
   _Ppc_Build( self, data );

}

void _Ppc_Auto_Thermal_Profile_Initialise( void* _self, void* data )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   Stg_Component_Initialise( self->manager, data, False );
   
   /* Initialize parent */
   _Ppc_Initialise( self, data );

}

void _Ppc_Auto_Thermal_Profile_Delete( void* _self )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_Auto_Thermal_Profile_Print( void* _self, Stream* stream )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_Auto_Thermal_Profile_Execute( void* _self, void* data )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_Auto_Thermal_Profile_Destroy( void* _self, void* data )
{
   Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_Auto_Thermal_Profile_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
	Ppc_Auto_Thermal_Profile* self = (Ppc_Auto_Thermal_Profile*) _self;
	Coord coord;
	double start_coord, end_coord, min_temp, max_temp, temp;
	double rho, cp, alpha, radiogenic;
	double k, c1, z, term1;
	unsigned int vertaxis;
	int err;
	
	err = PpcManager_Get( self->manager, lElement_I, particle, self->cpTag, &cp );
	assert(!err);
	err = PpcManager_Get( self->manager, lElement_I, particle, self->diffusivityTag, &alpha );
	assert(!err);
	err = PpcManager_Get( self->manager, lElement_I, particle, self->radiogenicTag, &radiogenic );
	assert(!err);
	err = PpcManager_Get( self->manager, lElement_I, particle, self->densityTag, &rho );
	assert(!err);
	
	start_coord = self->start_coord;
	end_coord = self->end_coord;
	min_temp = self->min_temp;
	max_temp = self->max_temp;
	vertaxis = self->axis;
	
    FeMesh_CoordLocalToGlobal( self->manager->mesh, lElement_I, particle->xi, coord );
    z = coord[vertaxis];
    
	if( z > start_coord )
		temp = min_temp;
	else if( z <= start_coord && z > end_coord ) {
		/* This is doing the steady-state conduction-advection geotherm */
		k = rho*cp*alpha;
		term1 = (radiogenic/(2*k))*(z*z);
		c1 = (max_temp-min_temp+term1)/z;
		temp = -term1 + (c1*z) + min_temp;	
	}
	else
		temp = max_temp;
	/*
    if( temp > max_temp )
        temp = max_temp;
    if( temp < min_temp )
        temp = min_temp;
    */    
    result[0] = temp;
    
	return 0;
}
