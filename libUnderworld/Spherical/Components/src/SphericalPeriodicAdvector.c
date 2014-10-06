#include <mpi.h>
#include <string.h>
#include <assert.h>

#include "Components.h"


char IJKTopology_DimNumToDimLetter[3] = {'I', 'J', 'K'};

/* Textual name of this class */
const Type SphericalPeriodicAdvector_Type = "SphericalPeriodicAdvector";

/* Constructors ------------------------------------------------------------------------------------------------*/
SphericalPeriodicAdvector* SphericalPeriodicAdvector_New( 
	Name						name,
	PICelleratorContext*	context,
	Mesh*						mesh, 
	Swarm*					swarm,
	Dictionary*				dictionary )
{
	SphericalPeriodicAdvector* self = _SphericalPeriodicAdvector_DefaultNew( name );

	self->isConstructed = True;
	_PeriodicBoundariesManager_Init( self, context, mesh, swarm, dictionary );
	_SphericalPeriodicAdvector_Init( self );

	return self;
}	

void* _SphericalPeriodicAdvector_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(SphericalPeriodicAdvector);
	Type                                                      type = SphericalPeriodicAdvector_Type;
	Stg_Class_DeleteFunction*                              _delete = _PeriodicBoundariesManager_Delete;
	Stg_Class_PrintFunction*                                _print = _PeriodicBoundariesManager_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SphericalPeriodicAdvector_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _SphericalPeriodicAdvector_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _SphericalPeriodicAdvector_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _PeriodicBoundariesManager_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _PeriodicBoundariesManager_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _PeriodicBoundariesManager_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;

	return (void*) _SphericalPeriodicAdvector_New(  SPHERICALPERIODICADVECTOR_PASSARGS  );
}

SphericalPeriodicAdvector* _SphericalPeriodicAdvector_New(  SPHERICALPERIODICADVECTOR_DEFARGS  ) {
	SphericalPeriodicAdvector* self;
	
	/* Allocate memory */
	self = (SphericalPeriodicAdvector*)_PeriodicBoundariesManager_New(  PERIODICBOUNDARIESMANAGER_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	self->_updateParticle = SphericalPeriodicAdvector_UpdateParticle;
	
	return self;
}

void _SphericalPeriodicAdvector_Init( SphericalPeriodicAdvector* self )
{
}

void _SphericalPeriodicAdvector_AssignFromXML( void* periodicBCsManager, Stg_ComponentFactory* cf, void* data ) {
	SphericalPeriodicAdvector*	self = (SphericalPeriodicAdvector*)periodicBCsManager;
	Dictionary*						dictionary = NULL;
	Mesh*								mesh = NULL;
	Swarm*							swarm = NULL;
	PICelleratorContext*			context;

	_PeriodicBoundariesManager_AssignFromXML( self, cf, data );

	_SphericalPeriodicAdvector_Init( self );
}

void _SphericalPeriodicAdvector_Delete( void* perBCsManager ) {
	SphericalPeriodicAdvector* self = (SphericalPeriodicAdvector*)perBCsManager;
	
	/* Stg_Class_Delete parent */
	_Stg_Component_Delete( self );
}

/* Virtual Functions -------------------------------------------------------------------------------------------------------------*/

void _SphericalPeriodicAdvector_Build( void* periodicBCsManager, void* data ) {	
	SphericalPeriodicAdvector* self = (SphericalPeriodicAdvector*)periodicBCsManager;
	Dictionary_Entry_Value*    periodicBCsList = NULL;

	Stg_Component_Build( self->swarm, data, False );
	Stg_Component_Build( self->mesh, data, False );
	self->size = 4;
	self->boundaries = Memory_Alloc_Array( PeriodicBoundary, self->size, "SphericalPeriodicAdvector->boundaries" );

	if ( self->dictionary ) {
		periodicBCsList = Dictionary_Get( self->dictionary, (Dictionary_Entry_Key)"PeriodicBoundaries" );
		
		/* Dictionary entry is optional - users may prefer to enter in code */
		if ( periodicBCsList ) {
			Index                   numPeriodicBCs = 0;
			Index                   periodicBC_I = 0;
			Dictionary_Entry_Value* periodicBC = NULL;
			char*                   perBCAxis = NULL;
			
			numPeriodicBCs = Dictionary_Entry_Value_GetCount( periodicBCsList );

			for ( periodicBC_I = 0; periodicBC_I < numPeriodicBCs; periodicBC_I++  ) {
				periodicBC = Dictionary_Entry_Value_GetElement( periodicBCsList, periodicBC_I );
				perBCAxis = Dictionary_Entry_Value_AsString( periodicBC );

				if ( 0 == strcmp( perBCAxis, "I_AXIS" ) ) {
					SphericalPeriodicAdvector_AddPeriodicBoundary( self, I_AXIS );
				}
				else if ( 0 == strcmp( perBCAxis, "J_AXIS" ) ) {
					SphericalPeriodicAdvector_AddPeriodicBoundary( self, J_AXIS );
				}
				else if ( 0 == strcmp( perBCAxis, "K_AXIS" ) ) {
					SphericalPeriodicAdvector_AddPeriodicBoundary( self, K_AXIS );
				}
			}
		}
	}
	/* Test if mesh is periodic */
	else if ( Stg_Class_IsInstance( self->mesh->generator, CartesianGenerator_Type ) ) {
		CartesianGenerator* cartesianGenerator = (CartesianGenerator*) self->mesh->generator;
		Dimension_Index dim_I;

		for ( dim_I = 0 ; dim_I < self->swarm->dim ; dim_I++ ) {
			/* Add boundaries straight from mesh generator */
			if ( cartesianGenerator->periodic[ dim_I ] ) 
				SphericalPeriodicAdvector_AddPeriodicBoundary( self, dim_I );
		}		
	}
}

/* Public Functions -------------------------------------------------------------------------------------------------------------*/

void SphericalPeriodicAdvector_AddPeriodicBoundary( void* periodicBCsManager, Axis axis ) {
	SphericalPeriodicAdvector*	self = (SphericalPeriodicAdvector*)periodicBCsManager;	
	PeriodicBoundary*				newPeriodicBoundary;
	double							min[3], max[3];

	Mesh_GetGlobalCoordRange( self->mesh, min, max );
	
	if ( self->count == self->size ) {
		self->size += self->delta;
		self->boundaries = Memory_Realloc_Array( self->boundaries, PeriodicBoundary, self->size );
	}

	newPeriodicBoundary = &self->boundaries[self->count];
	newPeriodicBoundary->axis = axis;
	newPeriodicBoundary->minWall = min[axis];
	newPeriodicBoundary->maxWall = max[axis];
	newPeriodicBoundary->particlesUpdatedMinEndCount = 0;	
	newPeriodicBoundary->particlesUpdatedMaxEndCount = 0;	
	self->count++;
}

void SphericalPeriodicAdvector_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I ) {
	Axis								boundaryAxis;	
	SphericalPeriodicAdvector*	self = (SphericalPeriodicAdvector*)periodicBCsManager;
	double							difference = 0.0;
	GlobalParticle*				particle = NULL;
	Index								perBoundary_I = 0;
	PeriodicBoundary*				perBoundary = NULL;

	Journal_DPrintfL( self->debug, 2, "In %s:\n", __func__ );
	Stream_Indent( self->debug );

	particle = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );

	Journal_DPrintfL( self->debug, 2, "Checking particle %d at (%.4g,%.4g,%.4g)\n", lParticle_I, particle->coord[0], particle->coord[1], particle->coord[2] );

	for ( perBoundary_I = 0; perBoundary_I < self->count; perBoundary_I++ ) {
		perBoundary = &self->boundaries[perBoundary_I];
		boundaryAxis = perBoundary->axis;

		Journal_DPrintfL( self->debug, 2, "Checking axis %d:\n", boundaryAxis );
		Stream_Indent( self->debug );

		if ( particle->coord[boundaryAxis] < perBoundary->minWall ) {
			Journal_DPrintfL( self->debug, 3, "coord is < min wall %.4f:\n", perBoundary->minWall );
			difference = perBoundary->minWall - particle->coord[boundaryAxis];
			particle->coord[boundaryAxis] = perBoundary->maxWall - difference;
			perBoundary->particlesUpdatedMinEndCount++;
			Journal_DPrintfL( self->debug, 3, "moving to (%.4f,%.4f,%.4f).\n", particle->coord[I_AXIS], particle->coord[J_AXIS], particle->coord[K_AXIS] );
		}
		else if ( particle->coord[perBoundary->axis] > perBoundary->maxWall ) {
			Journal_DPrintfL( self->debug, 3, "coord is > max wall %.4f:\n", perBoundary->maxWall );
			difference = particle->coord[boundaryAxis] - perBoundary->maxWall; 
			particle->coord[boundaryAxis] = perBoundary->minWall + difference;
			perBoundary->particlesUpdatedMaxEndCount++;
			Journal_DPrintfL( self->debug, 3, "moving to (%.4f,%.4f,%.4f).\n", particle->coord[I_AXIS], particle->coord[J_AXIS], particle->coord[K_AXIS] );
		}
		Stream_UnIndent( self->debug );
	}	

	Stream_UnIndent( self->debug );

	/* TODO: this is a bit of a hack to print this here using the lParticleI = swarm->total - 1, but its
	the only way I can see given this func is part of the SwarmAdvector intermediate. Should really be a 
	function on this class that updates all the particles. -- Main.PatrickSunter 15 May 2006 */
	if ( lParticle_I == (self->swarm->particleLocalCount-1) ) {
		PeriodicBoundary*	boundary = NULL;
		Index					perB_I;
	
		Journal_DPrintfL( self->debug, 1, "SphericalPeriodicAdvector total particles updated:\n" );
		Stream_Indent( self->debug );

		for ( perB_I = 0; perB_I < self->count; perB_I++ ) {
			boundary = &self->boundaries[perB_I];

			Journal_DPrintfL( self->debug, 1, "Periodic Boundary in %c Axis {%.2f,%.2f}: %d min end, %d max end\n",
				IJKTopology_DimNumToDimLetter[boundary->axis], boundary->minWall, boundary->maxWall,
				boundary->particlesUpdatedMinEndCount, boundary->particlesUpdatedMaxEndCount );
			/* Reset the counters for next time */
			boundary->particlesUpdatedMinEndCount = 0;	
			boundary->particlesUpdatedMaxEndCount = 0;	
		}
		Stream_UnIndent( self->debug );
	}
}	



