/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005-2010, Monash University 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%  Louis.Moresi - Louis.Moresi@monash.edu
*%
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "ViscoelasticRheology.h"
#include "ViscoelasticForceTerm.h"
#include "JaumannRotator.h"

#include <string.h>
#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type ViscoelasticForceTerm_Type = "ViscoelasticForceTerm";

void* _ViscoelasticForceTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(ViscoelasticForceTerm);
	Type                                                      type = ViscoelasticForceTerm_Type;
	Stg_Class_DeleteFunction*                              _delete = _ForceTerm_Delete;
	Stg_Class_PrintFunction*                                _print = _ForceTerm_Print;
	Stg_Class_CopyFunction*                                  _copy = _ForceTerm_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ViscoelasticForceTerm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ViscoelasticForceTerm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ForceTerm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ViscoelasticForceTerm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ForceTerm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ViscoelasticForceTerm_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _ViscoelasticForceTerm_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _ViscoelasticForceTerm_New(  VISCOELASTICFORCETERM_PASSARGS  );
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
ViscoelasticForceTerm* _ViscoelasticForceTerm_New(  VISCOELASTICFORCETERM_DEFARGS  ) 
{
	ViscoelasticForceTerm* self;

	assert( _sizeOfSelf >= sizeof(ViscoelasticForceTerm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (ViscoelasticForceTerm*) _ForceTerm_New(  FORCETERM_PASSARGS  );
	
	return self;
}

void _ViscoelasticForceTerm_Init(
	void*							forceTerm,
	FiniteElementContext*   context,
	ConstitutiveMatrix*     constitutiveMatrix,
	FeVariable*             strainRateField,
	MaterialPointsSwarm*    materialPointsSwarm,
	JaumannRotator*         jaumannRotator,
	Materials_Register*     materials_Register )
{
	ViscoelasticForceTerm*	self = (ViscoelasticForceTerm*)forceTerm;
	Viscoelastic_Particle*	particleExt;
	StandardParticle		particle;
	Dimension_Index			dim = materialPointsSwarm->dim;
	SwarmVariable*			swarmVariable;
	SwarmVariable*			stressSwarmVariable;
				
	/* Assign Pointers */
	self->dim = dim;
	self->constitutiveMatrix = constitutiveMatrix;

	self->strainRateField = strainRateField;
	self->materials_Register = materials_Register;
	self->materialPointsSwarm = materialPointsSwarm;
	self->jaumannRotator = jaumannRotator;
	
	/* Add function to AbstractContext_EP_UpdateClass entry point - which is run at the end of each step */
   /* FORBEC: You want to change the entry point this function _UpdateStress happens on */ 
	EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ), _ViscoelasticForceTerm_UpdateStress, self );

	/* Add particle extension */
	
	self->particleExtHandle  = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, (Name)ViscoelasticForceTerm_Type, sizeof(Viscoelastic_Particle)  );	
	
	/* Add SwarmVariables for plotting */
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &particle, self->particleExtHandle );
	
	swarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"PreviousStress", (ArithPointer) &particleExt->prevStress - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim),
		"Ps0","Ps1","Ps2","Ps3","Ps4","Ps5");
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );
	
	swarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"ElasticStress", (ArithPointer) &particleExt->elasticStress - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim),
		"Es0","Es1","Es2","Es3","Es4","Es5");
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );

	stressSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"Stress", (ArithPointer) &particleExt->totalStress - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		"St0","St1","St2","St3","St4","St5" );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)stressSwarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)stressSwarmVariable->variable );

	stressSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"StressRate", (ArithPointer) &particleExt->stressRate - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		"Sr0","Sr1","Sr2","Sr3","Sr4","Sr5" );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)stressSwarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)stressSwarmVariable->variable );
   	
	swarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"ParticleStrainRate", (ArithPointer) &particleExt->ParticleStrainRate - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		"SRate0","SRate1","SRate2","SRate3","SRate4","SRate5" );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );

   	swarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"ParticleViscousStrainRate", 
   		(ArithPointer) &particleExt->ParticleViscousStrainRate - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		"VSRate0","VSRate1","VSRate2","VSRate3","VSRate4","VSRate5" );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );

   	swarmVariable =  Swarm_NewScalarVariable( materialPointsSwarm, (Name)"ParticleViscousStrainRateInv", 
   		(ArithPointer) &particleExt->ParticleViscousStrainRateInv - (ArithPointer) &particle, 
   		Variable_DataType_Double  );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );

   	swarmVariable =  Swarm_NewScalarVariable( materialPointsSwarm, (Name)"ParticleOriginalViscosity", 
   		(ArithPointer) &particleExt->ParticleOriginalViscosity - (ArithPointer) &particle, 
   		Variable_DataType_Double  );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable );
   	LiveComponentRegister_IfRegThenAdd( (Stg_Component*)swarmVariable->variable );


	/* if using jaumannRotator, need to submit variable to be integrated in time */
	if ( jaumannRotator )
		jaumannRotator->variable = stressSwarmVariable->variable;
}

void _ViscoelasticForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ){
	ViscoelasticForceTerm*	self = (ViscoelasticForceTerm*)forceTerm;
	FeVariable *strainRateField;
	Materials_Register*		materials_Register;
	MaterialPointsSwarm*		materialPointsSwarm;
	JaumannRotator*			jaumannRotator;
	PICelleratorContext*		context;
	ConstitutiveMatrix*         constitutiveMatrix;
	
	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );
	
	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
   
	// TODO : KeyFallback deprecated soon

	strainRateField = Stg_ComponentFactory_ConstructByNameWithKeyFallback( cf, self->name, (Name)"StrainRateField", (Dictionary_Entry_Key)"StrainRateField", FeVariable, True, data  );
	jaumannRotator = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"JaumannRotator", JaumannRotator, False, data );
	constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ConstitutiveMatrix", ConstitutiveMatrix, True, data );
	context = (PICelleratorContext* )self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register;
	assert( materials_Register );
	
	_ViscoelasticForceTerm_Init(
		self ,
		(FiniteElementContext*) context,
		constitutiveMatrix,
		strainRateField, 
		materialPointsSwarm,
		jaumannRotator,
		materials_Register );
}


void _ViscoelasticForceTerm_Initialise( void* forceTerm, void* data ) {
   
   ViscoelasticForceTerm* self = (ViscoelasticForceTerm*) forceTerm;
   MaterialPointsSwarm *materialSwarm = self->materialPointsSwarm;
   MaterialPoint *materialPoint;
   Material *material = NULL;
   Viscoelastic_Particle *pve_ext = NULL;
   int particle_I, dofCount, handle_pev_ext;
   int particleLocalCount = materialSwarm->particleLocalCount;
   int somethingIsCompressible = 0;

   handle_pev_ext = self->particleExtHandle;
   dofCount = StGermain_nSymmetricTensorVectorComponents(self->dim);
	
	_ForceTerm_Initialise( self, data );
	
	if ( self->context->loadSwarmsFromCheckpoint == False ) {
		
	   if( self->materialPointsSwarm )				Stg_Component_Build( self->materialPointsSwarm, data, False );
	   if( self->constitutiveMatrix ) 				Stg_Component_Build( self->constitutiveMatrix, data, False ); 
	   if( self->strainRateField )    				Stg_Component_Build( self->strainRateField, data, False );
	   if( self->ForceStress )        				Stg_Component_Build( self->ForceStress, data, False );
	   if( self->Stress )							Stg_Component_Build( self->Stress, data, False );
	   if( self->ParticleStrainRate ) 				Stg_Component_Build( self->ParticleStrainRate, data, False );
	   if( self->ParticleViscousStrainRate )   		Stg_Component_Build( self->ParticleViscousStrainRate, data, False );
	   if( self->ParticleViscousStrainRateInv )		Stg_Component_Build( self->ParticleViscousStrainRateInv, data, False );
	   if( self->ParticleOriginalViscosity )		Stg_Component_Build( self->ParticleOriginalViscosity, data, False );
	   if( self->jaumannRotator )					Stg_Component_Build( self->jaumannRotator, data, False );

	   for ( particle_I = 0 ; particle_I < particleLocalCount ; particle_I++ ) {
			materialPoint    = (MaterialPoint*)Swarm_ParticleAt( materialSwarm, particle_I );
			material = Materials_Register_GetByIndex( materialSwarm->materials_Register, materialPoint->materialIndex );

	      /* get viscoelastic particle extension */
	      pve_ext = ExtensionManager_Get( materialSwarm->particleExtensionMgr, materialPoint, handle_pev_ext );
	      assert( pve_ext != NULL );

	      /* zero the initial prevStress pointer */ 
	      memset( pve_ext->prevStress, 0, dofCount*sizeof(double) );

	      /* set the ParticleOriginalViscosity to zero (so we can test it) */
	      pve_ext->ParticleOriginalViscosity = 0.0;

	   }
	}	
}

void _ViscoelasticForceTerm_Destroy( void* _self, void* data ) {
	ViscoelasticForceTerm* self = (ViscoelasticForceTerm*) _self;

   if( self->materialPointsSwarm )				Stg_Component_Destroy( self->materialPointsSwarm, data, False );
   if( self->constitutiveMatrix )				Stg_Component_Destroy( self->constitutiveMatrix, data, False ); 
   if( self->strainRateField )					Stg_Component_Destroy( self->strainRateField, data, False );
   if( self->ForceStress )						Stg_Component_Destroy( self->ForceStress, data, False );
   if( self->Stress )							Stg_Component_Destroy( self->Stress, data, False );
   if( self->ParticleStrainRate )				Stg_Component_Destroy( self->ParticleStrainRate, data, False );
   if( self->ParticleViscousStrainRate )		Stg_Component_Destroy( self->ParticleViscousStrainRate, data, False );
   if( self->ParticleViscousStrainRateInv )		Stg_Component_Destroy( self->ParticleViscousStrainRateInv, data, False );
   if( self->ParticleOriginalViscosity )		Stg_Component_Destroy( self->ParticleOriginalViscosity, data, False );
   if( self->jaumannRotator ) 					Stg_Component_Destroy( self->jaumannRotator, data, False );

   _ForceTerm_Destroy( self, data ); 
}

void _ViscoelasticForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
    ViscoelasticForceTerm*    			self                = (ViscoelasticForceTerm*) forceTerm;
    IntegrationPoint*         			integrationPoint;
    MaterialPoint*            			materialPoint;
    double*                   			xi;
    Particle_InCellIndex      			cParticle_I;
    Particle_InCellIndex      			cellParticleCount;
    Element_NodeIndex         			elementNodeCount;
    Node_ElementLocalIndex    			node_I;
    ElementType*              			elementType;
    Dof_Index                 			dofsPerNode;
    Cell_Index                			cell_I;
    double**                  			GNx;
    double                    			detJac;
    double                    			weight;
    double                    			factor;
    double                    			mu;
    Viscoelastic_Particle*    			pve_ext;
    Material*                 			material;
    Dimension_Index           			dim                 = self->dim;
    double                    			dt_e;
    double                    			eta_eff;
    double                    			eta;
    IntegrationPointsSwarm*   			integrationSwarm    = (IntegrationPointsSwarm*)self->integrationSwarm;
    MaterialPointsSwarm*      			materialPointsSwarm = self->materialPointsSwarm;
    ViscoelasticRheology*     			viscoelasticRheology;
    FeMesh*                   			mesh                = forceVector->feVariable->feMesh;
    SymmetricTensor strainRate;
    int dofCount, dof_i;
	SymmetricTensor totalStress;
    double *tau_prev;

    /* handle for the stored constitutive matrix on material particles */
    int handle_storedmatrix;

    /* handle for the particle viscoelastic extension */
    int handle_pev_ext = self->particleExtHandle;

    /* get the dofCount for a 2nd order symmetric tensor */
    dofCount = StGermain_nSymmetricTensorVectorComponents( dim ); 

    /* if the constitutive matrix doesn't store particle const. matrices
     * set the handle to -1 */

    if( False == self->constitutiveMatrix->storeConstitutiveMatrix )
		handle_storedmatrix = -1;
    else
		handle_storedmatrix = ConstitutiveMatrix_GetParticleConstExtHandle( self->constitutiveMatrix );
	
    /* Set the element type */
    elementType      = FeMesh_GetElementType( mesh, lElement_I );
    elementNodeCount = elementType->nodeCount;	

    /* assumes constant number of dofs per element */
    dofsPerNode = dim;
		
    cell_I = CellLayout_MapElementIdToCellId( integrationSwarm->cellLayout, lElement_I );
    cellParticleCount = integrationSwarm->cellParticleCountTbl[ cell_I ];

    GNx = Memory_Alloc_2DArray( double, dim, elementNodeCount, (Name)"GNx" );
    double *cm = NULL;

    for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		integrationPoint = (IntegrationPoint* ) Swarm_ParticleInCellAt( integrationSwarm, cell_I, cParticle_I );
		material = IntegrationPointsSwarm_GetMaterialOn( integrationSwarm, integrationPoint );
		
		/* check if the material data structure is of a RheologyMaterial_Type */
		if( !Stg_Class_IsInstance ( material, RheologyMaterial_Type ) )
		    continue;

		/* Check if this particle is viscoelastic */
		viscoelasticRheology = (ViscoelasticRheology*) RheologyMaterial_GetRheologyByType( material, ViscoelasticRheology_Type );
	
		if ( ! viscoelasticRheology ) 
		    continue;

		xi      = integrationPoint->xi;
		weight  = integrationPoint->weight;
		mu      = viscoelasticRheology->mu;
		dt_e    = viscoelasticRheology->elasticTimeStep;

		ElementType_ShapeFunctionsGlobalDerivs( elementType, mesh, lElement_I, xi, dim, &detJac, GNx );

		/* Ok, we break any fancy mapper stuff here and assume the materialParticles correspond exactly to the Integration ones */
		materialPoint = OneToOneMapper_GetMaterialPoint( integrationSwarm->mapper, integrationPoint, &materialPointsSwarm );

		/* get the particle viscoelastic extension */
		pve_ext = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, handle_pev_ext );

		/* get the isotropic viscosity off the material point.
		 * Either we need to assemble the constitutive matrix again or we call the storedmatrix extension */

		if( handle_storedmatrix == -1 ) {
		    ConstitutiveMatrix_Assemble( self->constitutiveMatrix, lElement_I, cParticle_I, integrationPoint ); 
		    eta_eff = ConstitutiveMatrix_GetIsotropicViscosity( self->constitutiveMatrix );
		} 
		else {
		    cm = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, handle_storedmatrix );
		    assert( cm != NULL );
		    if (dim == 2 ) 
				eta_eff = cm[8];
		    else 
				eta_eff = cm[21];
		}

		/* This should have been stored when the viscoelasticRheology was called */
		eta = eta_eff * mu * dt_e / ( mu * dt_e - eta_eff);
	
		tau_prev = pve_ext->prevStress;

		/* For analysis purposes, we update these particle extensions */

		FeVariable_InterpolateWithinElement( self->strainRateField, materialPoint->owningCell, xi, strainRate );

		for( dof_i = 0 ; dof_i < dofCount ; dof_i++ ) {
            /* save "elastic stress", for this timestep, on particle extension */
            pve_ext->elasticStress[dof_i] = tau_prev[dof_i];

			totalStress[dof_i] = 2*eta_eff*strainRate[dof_i] + eta_eff/(mu*dt_e)*tau_prev[dof_i];
            pve_ext->stressRate[dof_i] = ( totalStress[dof_i] - pve_ext->totalStress[dof_i] ) / dt_e; //forceVector->context->dt ;
			pve_ext->ParticleStrainRate[dof_i] = strainRate[dof_i];
            /* save tau_de on the particle extension - for viz of total stress for this timestep */
            pve_ext->totalStress[dof_i] = totalStress[dof_i];     
            pve_ext->ParticleViscousStrainRate[dof_i] = 0.5 * totalStress[dof_i] / eta ;    
		}

		pve_ext->ParticleViscousStrainRateInv = 0.5 * SymmetricTensor_2ndInvariant( totalStress, self->dim ) / eta;
  
		/* integration weight multiplied by eta_eff/(mu*dt_e) */
		assert( fabs(mu*dt_e) > 1e-15 );
		factor = detJac * integrationPoint->weight * eta_eff / (mu*dt_e); 

		/* construct force */

		/* 	Note: 2D contains terms related to compressibility, 3D does not. This needs to be fixed
	   		and I am not sure if ViscoelasticCompressibleForceTerm.c is then not needed - LM */

		if( dim == 2 ) {
		    for( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) {
				elForceVec[node_I*dofsPerNode + I_AXIS ] -= factor * 
				    (GNx[0][node_I]*tau_prev[0] + GNx[1][node_I]*tau_prev[2] ) ;
				elForceVec[node_I*dofsPerNode + J_AXIS ] -= factor * 
				    (GNx[1][node_I]*tau_prev[1] + GNx[0][node_I]*tau_prev[2] ) ;
		    }
		} 
		else {
		    for( node_I = 0 ; node_I < elementNodeCount ; node_I++ ) {
			elForceVec[node_I*dofsPerNode + I_AXIS ] -= factor * 
			    (GNx[0][node_I]*tau_prev[0] + GNx[1][node_I]*tau_prev[3] + GNx[2][node_I]*tau_prev[4]);

			elForceVec[node_I*dofsPerNode + J_AXIS ] -= factor * 
			    (GNx[1][node_I]*tau_prev[1] + GNx[0][node_I]*tau_prev[3] + GNx[2][node_I]*tau_prev[5] );

			elForceVec[node_I*dofsPerNode + K_AXIS ] -= factor * 
			    (GNx[2][node_I]*tau_prev[2] + GNx[0][node_I]*tau_prev[4] + GNx[1][node_I]*tau_prev[5] );
		    }
		}
    }

    Memory_Free(GNx); 
}

void _ViscoelasticForceTerm_UpdateStress( void* forceTerm, void* data ) {
	ViscoelasticForceTerm*            self               = (ViscoelasticForceTerm*) forceTerm;
	PICelleratorContext*                 context            = (PICelleratorContext*) data;
	double                               dt_e;
	double                               dt                 = context->dt;
	double                               phi;
	MaterialPoint*                       materialPoint;
	Viscoelastic_Particle*            	 pve_ext;
	MaterialPointsSwarm*                 materialPointsSwarm = (MaterialPointsSwarm*)self->materialPointsSwarm;
	Dimension_Index                      dim                = self->dim;
	Material*                            material;
	Dof_Index                            dofCount, dof_i;
	Particle_Index                       particle_I         = 0;
	double                               eta_eff;
	ViscoelasticRheology*                viscoelasticRheology;
	double                               mu;
	Coord                                xi;
	
	
	int    localParticleCount=0;
	int    globalParticleCount=0;
	double particleStoredStressInvariant;
	double localMeanStoredStressInvariant=0.0;
	double globalMeanStoredStressInvariant=0.0;	


   SymmetricTensor strainRate, tau_de, tau_dt;
   double *cm = NULL;
   int ii, cell_i, particleCellCount, lCellCount;

   /* handle for the stored constitutive matrix on material particles */
   int handle_storedmatrix;

   /* handle for the particle viscoelastic extension */
   int handle_pev_ext = self->particleExtHandle;
	
   /* if the constitutive matrix doesn't store particle const. matrices
    * set the handle to -1 */

   if( False == self->constitutiveMatrix->storeConstitutiveMatrix )
      handle_storedmatrix = -1;
   else
      handle_storedmatrix = ConstitutiveMatrix_GetParticleConstExtHandle( self->constitutiveMatrix );

   /* check if this swarm layout is compatible with this algorithm which
    * assumes cells and elements map to one another */
   assert( Stg_CheckType( materialPointsSwarm->cellLayout, ElementCellLayout ) );

	Journal_DPrintf( self->debug, "In %s for %s '%s'\n", __func__, self->type, self->name );
	
   /* get the dofCount for a 2nd order symmetric tensor */
	dofCount = StGermain_nSymmetricTensorVectorComponents(dim);
   lCellCount = materialPointsSwarm->cellLocalCount;
	
   for( cell_i = 0 ; cell_i < lCellCount ; cell_i++ ) {
      particleCellCount = materialPointsSwarm->cellParticleCountTbl[cell_i];

      for ( particle_I = 0 ; particle_I < particleCellCount ; particle_I++ ) {
         materialPoint = (MaterialPoint*)Swarm_ParticleInCellAt( materialPointsSwarm, cell_i, particle_I );
         material = Materials_Register_GetByIndex( materialPointsSwarm->materials_Register, materialPoint->materialIndex );

         /* check if the material data structure is of a RheologyMaterial_Type */
         if( !Stg_Class_IsInstance ( material, RheologyMaterial_Type ) )
            continue;

         /* Check if this particle is viscoelastic */
         viscoelasticRheology = (ViscoelasticRheology*)  RheologyMaterial_GetRheologyByType( material, ViscoelasticRheology_Type );
         if ( ! viscoelasticRheology ) 
               continue;
               
        dt_e = viscoelasticRheology->elasticTimeStep; 
		assert( fabs(dt_e) > 1e-15 );
        phi = dt / dt_e;
        mu = viscoelasticRheology->mu;
         
         /* get viscoelastic particle extension */
         pve_ext = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, handle_pev_ext );

         /* for interpolation of mesh quantity we need the local element coords */
         FeMesh_CoordGlobalToLocal( self->strainRateField->feMesh, materialPoint->owningCell, materialPoint->coord, xi );

         /* find strain rate at particle location */
         FeVariable_InterpolateWithinElement( self->strainRateField, materialPoint->owningCell, xi, strainRate );

         /* construct total deviatoric stress (=viscous stress + elastic stress) and save in particleExt->prevStress */

         /* get the isotropic viscosity off the material point.
          * Either we need to assemble the constitutive matrix again or we call the storedmatrix extension */

         if( handle_storedmatrix == -1 ) {
            ConstitutiveMatrix_AssembleMaterialPoint( self->constitutiveMatrix, materialPoint->owningCell, materialPointsSwarm, materialPointsSwarm->cellParticleTbl[cell_i][particle_I] ); 
            eta_eff = ConstitutiveMatrix_GetIsotropicViscosity( self->constitutiveMatrix );
         } 
		 else {  /* This does not get correctly stored during a checkpoint step */
            cm = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, handle_storedmatrix );
            assert( cm != NULL );
            if (dim == 2 ) 
				eta_eff = cm[8];
            else 
				eta_eff = cm[21];
         }

         assert( fabs(mu) > 1e-15  );
         /* calculate deviatoric de stress */
         for( dof_i = 0 ; dof_i < dofCount ; dof_i++ ) {
            tau_de[dof_i] = 
               2*eta_eff*strainRate[dof_i] + 
               eta_eff/(mu*dt_e)*pve_ext->prevStress[dof_i];
         }

         
         /* smooth de stress with pve->prevStress */
         for( ii = 0 ; ii < dofCount; ii++ ) {
            tau_dt[ii] = phi*tau_de[ii] + (1-phi)*pve_ext->prevStress[ii]; 
			// fprintf(stderr,"stored = %d, eta_eff = %g; tau[%d] = %g*%g + %g*%g\n", handle_storedmatrix, eta_eff,ii,phi,tau_de[ii],(1-phi),pve_ext->prevStress[ii]);
         }

         /* copy smoothed stress into the pve->prevStress buffer */
         memcpy( pve_ext->prevStress, tau_dt, dofCount*sizeof(double) );

         if (dim == 2 ) 
			particleStoredStressInvariant = 
				tau_dt[0] * tau_dt[0] + 
				tau_dt[1] * tau_dt[1] + 
				2 * tau_dt[2] * tau_dt[2];
		 else
			particleStoredStressInvariant = 
				tau_dt[0] * tau_dt[0] + 
				tau_dt[1] * tau_dt[1] + 
				tau_dt[2] * tau_dt[2] + 
				2 * tau_dt[3] * tau_dt[3] +
				2 * tau_dt[4] * tau_dt[4] +
				2 * tau_dt[5] * tau_dt[5];
						
		particleStoredStressInvariant = sqrt(particleStoredStressInvariant); 		 
		localMeanStoredStressInvariant += particleStoredStressInvariant;
		localParticleCount += 1;

      }
 	}

	/* These are the averaged values over particles in the domain which are viscoelastic */
	(void)MPI_Allreduce( &localMeanStoredStressInvariant,&globalMeanStoredStressInvariant, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	(void)MPI_Allreduce( &localParticleCount,&globalParticleCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	
	globalMeanStoredStressInvariant /= globalParticleCount;
	
	/* Note: there appears to be a bug with the stored matrix version of this code (eta_eff = 0 after restart)*/
	
	PetscPrintf(MPI_COMM_WORLD, "Mean Stored Stress Invariant = %g (on %d particles)\n", globalMeanStoredStressInvariant, globalParticleCount);
	//PetscPrintf(MPI_COMM_WORLD, "Effective viscosity = %g; dT = %g; dTe = %g; Mu = %g; phi = %g\n",eta_eff,dt,dt_e, mu, phi );

}



