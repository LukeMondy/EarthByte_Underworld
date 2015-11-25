
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Viscoelastic/Base/Base.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>


#include <assert.h>
#include <string.h>

Particle_Index            closestParticle_I;

const Type Viscoelastic_ViscoelasticCantileverBeam_Type = "Viscoelastic_ViscoelasticCantileverBeam";

void TurnOffGravity( UnderworldContext* context ) {
	BuoyancyForceTerm* buoyancy;
   Stokes_SLE* stokesSLE = (Stokes_SLE* )LiveComponentRegister_Get( context->CF->LCRegister, (Name)"stokesEqn" );

	buoyancy = (BuoyancyForceTerm*) Stg_ObjectList_Get( stokesSLE->fForceVec->forceTermList, (Name)"buoyancyForceTerm" );
	assert( buoyancy );
	buoyancy->gravity = 0.0;
}

void CalcDeflectionAndCheckGravity( UnderworldContext* context ) {

	IntegrationPointsSwarm*    swarm              = (IntegrationPointsSwarm*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"picIntegrationPoints" );
	GlobalParticle*            particle; 

	MaterialPointsSwarm**      materialSwarms;
	MaterialPointsSwarm*       materialSwarm;
	Index                      swarmLength;

	double        defMax            = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"defMax", 0.35  );
	double        deflection        = 0.0;
	Coord         coord;
	double        distance;
	static int    visits            = 0;
	static double initial_position  = 0.0;

	materialSwarms = IntegrationPointMapper_GetMaterialPointsSwarms( swarm->mapper, &swarmLength );
	/* assume first swarm for now!!! */
	materialSwarm = materialSwarms[0];

	if (visits == 0) {
		coord[0] = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"x_mid", 0.75  );
		coord[1] = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"y_mid", 0.75  );
		coord[2] = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"z_mid", 0.0  );
		closestParticle_I = Swarm_FindClosestParticle( materialSwarm, 2, coord, &distance );
	}

	particle = (GlobalParticle*)Swarm_ParticleAt( materialSwarm, closestParticle_I );

	if (visits==0)
		initial_position = particle->coord[1];
	
	visits=1;

	deflection = initial_position - particle->coord[1];
	
	if ( deflection >= defMax ) {
		TurnOffGravity( context );
	}

	StgFEM_FrequentOutput_PrintValue( context, deflection );
}

void _Viscoelastic_ViscoelasticCantileverBeam_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
	
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, CalcDeflectionAndCheckGravity );

	StgFEM_FrequentOutput_PrintString( context, "Deflection" );
}

/* This function will provide StGermain the abilty to instantiate (create) this codelet on demand. */
void* _Viscoelastic_ViscoelasticCantileverBeam_DefaultNew( Name name ) {
	return Codelet_New(
			Viscoelastic_ViscoelasticCantileverBeam_Type,
			_Viscoelastic_ViscoelasticCantileverBeam_DefaultNew,
			_Viscoelastic_ViscoelasticCantileverBeam_AssignFromXML,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

/* This function is automatically run by StGermain when this plugin is loaded. The name must be "<plugin-name>_Register". */
Index Viscoelastic_ViscoelasticCantileverBeam_Register( PluginsManager* pluginsManager ) {
	/* A plugin is only properly registered once it returns the handle provided when submitting a codelet to StGermain. */
	return PluginsManager_Submit( pluginsManager, Viscoelastic_ViscoelasticCantileverBeam_Type, (Name)"0", _Viscoelastic_ViscoelasticCantileverBeam_DefaultNew  );
}


