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

/********************************************************************************************************
 For now this depends on ViscoelasticRheology being setup and inserted in a rheology list just before it,
 to set correct form of effective viscosity for visco-elasticity.

 It turns out that the Jacobian terms look identical to the usual DruckerPrager jacobian except that 
 the strain-rate is modified with the extra stored previous stress terms.

 ********************************************************************************************************/
#include <math.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
//#include "RheologyClass.h"
//#include "StrainWeakening.h"
//#include "YieldRheology.h"
#include "ViscoelasticForceTerm.h"
#include "ViscoelasticRheology.h"
#include "DruckerPragerXXVE.h"
//#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type DruckerPragerXXVE_Type = "DruckerPragerXXVE";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
DruckerPragerXXVE* _DruckerPragerXXVE_New(  DRUCKERPRAGERXXVE_DEFARGS  ) 
{
      DruckerPragerXXVE*					self;

      /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
      assert( _sizeOfSelf >= sizeof(DruckerPragerXXVE) );
      self = (DruckerPragerXXVE*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
      /* Function pointers for this class that are not on the parent class should be set here */
	
      return self;
}

void _DruckerPragerXXVE_Init( 
      DruckerPragerXXVE*               self,
      MaterialPointsSwarm*     materialPointsSwarm,
      IntegrationPointsSwarm*  integrationSwarm,
      FiniteElementContext*    context,
      FeVariable*              strainRateField,
      SwarmVariable*           swarmStrainRate,
      FeVariable*              pressureField, 
      SwarmVariable*           swarmPressure, 
      double                   cohesion,
      double                   frictionCoeff,
      double                   minYieldCriterion,
      double                   minViscosity,
      double                   alpha,
      double                   beta,
      double                   kappa,
      double                   alphaB,
      double                   kappaB,
      double                   softeningStrain,
      double                   initialDamageFraction,
      double                   initialDamageWavenumber,
      double                   initialDamageWavenumberSinI,
      double                   initialDamageWavenumberCosI,
      double                   initialDamageWavenumberSinK,
      double                   initialDamageWavenumberCosK,
      double                   initialDamageFactor,
      long int                 randomSeed,
      Stg_Shape*               initialStrainShape  )
{
      DruckerPragerXXVE_ParticleExt*  particleExt;
      StandardParticle materialPoint;

      self->strainRateField        = strainRateField;
      self->swarmStrainRate        = swarmStrainRate;
      self->cohesion               = cohesion;

      self->pressureField          = pressureField;
      self->swarmPressure          = swarmPressure;
      self->frictionCoeff          = frictionCoeff;
      
      self->materialPointsSwarm    = materialPointsSwarm;
      self->integrationSwarm       = integrationSwarm;
      self->context                = (PICelleratorContext*)context;
      self->minYieldCriterion      = minYieldCriterion;  /* truncate yield criterion  */
      self->minViscosity           = minViscosity;
      self->alpha                  = alpha; /* multiplier for strain in cohesion function = 1.0/softeningStrain */
      self->beta                   = beta;  /* healing factor: set to 0 for no healing */
      self->kappa                  = kappa; /* ratio of lower cohesion bound */
      self->alphaB                 = alphaB; /* multiplier for strain in frictionCoeff function = alpha now */
      self->kappaB                 = kappaB; /* ratio of lower frictionCoeff bound */

      self->softeningStrain             = softeningStrain;
      self->initialDamageFraction       = initialDamageFraction;
      self->initialDamageWavenumber     = initialDamageWavenumber;
      self->initialDamageWavenumberSinI = initialDamageWavenumberSinI ;
      self->initialDamageWavenumberCosI = initialDamageWavenumberCosI;
      self->initialDamageWavenumberSinK = initialDamageWavenumberSinK;
      self->initialDamageWavenumberCosK = initialDamageWavenumberCosK ;
      self->initialDamageFactor         = initialDamageFactor;
      self->randomSeed                  = randomSeed;
      /* calculate machine epsilon - using this for a tolerance for when strain-rate invariant approaches zero */
      self->eps = 1.0;
      self->epsI = 0;

      while( self->eps > 0.0){
	       self->eps = self->eps/2.0;
	       self->epsI++;
      }

      self->particleExtHandle = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, (Name)DruckerPragerXXVE_Type, sizeof(DruckerPragerXXVE_ParticleExt) );

      EP_PrependClassHook( Context_GetEntryPoint( context, AbstractContext_EP_PostSolvePreUpdateClass ), DruckerPragerXXVE_IntegrateStrain, self );

      /* get the extension */
      particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &materialPoint, self->particleExtHandle );
      /* Setup a plasticStrain Variable for Visualisation */
      self->plasticStrain = Swarm_NewScalarVariable( materialPointsSwarm, (Name)"DruckerPragerXXVEPlasticStrain", (ArithPointer) &particleExt->pstrain - (ArithPointer) &materialPoint, Variable_DataType_Double  );
      /*  */
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->plasticStrain->variable );
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->plasticStrain );

      /* stick an update variable function on the dumpclass EP */
      EP_PrependClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_DumpClass ), _DruckerPragerXXVE_UpdateDrawingVariables, self );
}

void _DruckerPragerXXVE_UpdateDrawingVariables( void* rheology ){
    DruckerPragerXXVE* self = (DruckerPragerXXVE*) rheology;
    MaterialPointsSwarm*     materialPointsSwarm = self->materialPointsSwarm;
    IntegrationPointsSwarm*  integrationSwarm    = self->integrationSwarm;
    DruckerPragerXXVE_ParticleExt*  particleExt;
    int elementIndex, particleIndex;

    Cell_LocalIndex  nElements = materialPointsSwarm->cellLocalCount;
    
    /* Update swarm variables */
    Variable_Update( self->plasticStrain->variable );
    
    for( elementIndex = 0; elementIndex < nElements; elementIndex++ ) {
	     int cellIndex  = CellLayout_MapElementIdToCellId( materialPointsSwarm->cellLayout, elementIndex );
	     int nParticles = materialPointsSwarm->cellParticleCountTbl[cellIndex];
	     for( particleIndex = 0; particleIndex < nParticles; particleIndex++ ) {

	       IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleInCellAt( integrationSwarm, cellIndex, particleIndex );
	       MaterialPoint*    matParticle = OneToOneMapper_GetMaterialPoint( integrationSwarm->mapper, integrationPoint, &(materialPointsSwarm) );
	       particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, matParticle, self->particleExtHandle );
	       Variable_SetValueFloat( self->plasticStrain->variable, particleIndex, particleExt->pstrain );
    	}
    }
}
static double integrateStrainOnParticle(){

    return 0.0;
}
void DruckerPragerXXVE_IntegrateStrain( void* rheology ){
      DruckerPragerXXVE* self = (DruckerPragerXXVE*) rheology;
      MaterialPointsSwarm*     materialPointsSwarm = self->materialPointsSwarm;
      IntegrationPointsSwarm*  integrationSwarm    = self->integrationSwarm;
      PICelleratorContext*     context             = self->context;
      double                   cohesion            = self->cohesion;
      double                   frictionCoeff       = self->frictionCoeff;
      double                   eta0                = self->eta0;  // should also have min Yield in here
      double                   alpha               = self->alpha;
      double                   beta                = self->beta;
      double                   kappa               = self->kappa;
      double                   alphaB              = self->alphaB;
      double                   kappaB              = self->kappaB;
      double                   dt                  = context->dt;
      double                   etavm, strainInv, c0,c1, yieldCriterion, yieldIndicator, b0,b1, pressure;
      Particle_Index           particleIndex, elementIndex;
      SymmetricTensor          strainRate;
      DruckerPragerXXVE_ParticleExt*  particleExt;
      
      if( !(self->strainRateField) ){
	    printf("Cannot find strainRateField in %s\n",__func__);
	    exit(1);      }
      if( !(self->pressureField) ){
	    printf("Cannot find pressureField in %s\n",__func__);
	    exit(1);      }
      //fprintf(stderr,"Integrate Strain\n");

      /* integrate the pstrain on the particles in time: heal if plastic SRI less than viscous SRI else pstrain increases (SRI= Strain-Rate second Invariant) */
      Cell_LocalIndex  nElements = materialPointsSwarm->cellLocalCount;
      for( elementIndex = 0; elementIndex < nElements; elementIndex++ ) {
    	    int cellIndex  = CellLayout_MapElementIdToCellId( materialPointsSwarm->cellLayout, elementIndex );
    	    int nParticles = materialPointsSwarm->cellParticleCountTbl[cellIndex];
      	    for( particleIndex = 0; particleIndex < nParticles; particleIndex++ ) {
          		  IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleInCellAt( integrationSwarm, cellIndex, particleIndex );
          		  MaterialPoint*    matParticle = OneToOneMapper_GetMaterialPoint( integrationSwarm->mapper, integrationPoint, &(materialPointsSwarm) );
          		  /* Get Parameters From Material Swarm Extension */
          		  particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, matParticle, self->particleExtHandle );
      		 
                FeVariable_InterpolateWithinElement( self->strainRateField, cellIndex, integrationPoint->xi, strainRate );


  		  /*****************************************************************************************************************/

  		  /* apply viscoelastic correction */

  		  if ( self->viscoelastic_corr ){
            /* e -> e' = e + 1/(2*mu*dte)*tau_0; tau_0 = previous stress */

  		      /* get viscoelastic particle extension */
  		      Viscoelastic_Particle* veExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, matParticle, self->veForceTerm->particleExtHandle );

  		      /* get viscoelastic rheology */
  		      ViscoelasticRheology* veRheology = self->veRheology;
  		      
  		      double mu = veRheology->mu;
  		      double dt_e = veRheology->elasticTimeStep;
  		      double factor = 0.5 / (mu*dt_e);
  		      int dim=context->dim;
  		      int i;
  		      for(i=0;i< dim*(dim+1)/2; i++){
  			         strainRate[i] = strainRate[i] + factor * veExt->prevStress[i]; 
            }
  		  }

  		  /*****************************************************************************************************************/
  		  
        strainInv = SymmetricTensor_2ndInvariant( strainRate, context->dim );/* is sqrt(I2) where I2 is strain-rate second invariant E_ij*E_ij/2 where E_ij is strain-rate tensor*/ 
  		  FeVariable_InterpolateWithinElement( self->pressureField, cellIndex, integrationPoint->xi, &pressure );
  		  
  		  double dg = 0.0;
  		  double Idt,gamma;
  		  gamma = particleExt->pstrain;
  		  c0 = cohesion;
  		  c1 = kappa*c0;
  		  b0 = frictionCoeff;
  		  b1 = kappaB*b0;
  		  Idt = strainInv*dt;
  		  
  		  yieldCriterion = (c0-c1)*exp(-(gamma)*alpha) + c1  + pressure*( (b0-b1)*exp(-(gamma)*alphaB) + b1) ;
  		  yieldIndicator = 2*eta0*strainInv;
  		  
  		  if(yieldIndicator > yieldCriterion){
  		      /* smooth monotonically decreasing curve for cohesion function */
  		      etavm = yieldCriterion/(2.0*strainInv);
  		      dg = (1.0 - (1.0+beta)*etavm/eta0)*strainInv*dt;
  		      particleExt->pstrain += dg;/* Update New pstrain for next time-step with Euler method */
  		      if(particleExt->pstrain < 0.0) 
  				particleExt->pstrain = 0.0; 
  		  }
  	    }
      }
      
}
void DruckerPragerXXVE_UpdateDt( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context ){
    //UnderworldContext*  context = (UnderworldContext*)_context;
    //FiniteElementContext_CalcNewDt( context );
}

void* _DruckerPragerXXVE_DefaultNew( Name name ) {
      /* Variables set in this function */
      SizeT                                                     _sizeOfSelf = sizeof(DruckerPragerXXVE);
      Type                                                             type = DruckerPragerXXVE_Type;
      Stg_Class_DeleteFunction*                                     _delete = _YieldRheology_Delete;
      Stg_Class_PrintFunction*                                       _print = _YieldRheology_Print;
      Stg_Class_CopyFunction*                                         _copy = _YieldRheology_Copy;
      Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _DruckerPragerXXVE_DefaultNew;
      Stg_Component_ConstructFunction*                           _construct = _DruckerPragerXXVE_AssignFromXML;
      Stg_Component_BuildFunction*                                   _build = _DruckerPragerXXVE_Build;
      Stg_Component_InitialiseFunction*                         _initialise = _DruckerPragerXXVE_Initialise;
      Stg_Component_ExecuteFunction*                               _execute = _YieldRheology_Execute;
      Stg_Component_DestroyFunction*                               _destroy = _DruckerPragerXXVE_Destroy;
      Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _DruckerPragerXXVE_ModifyConstitutiveMatrix;

      /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
      AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

      return (void*) _DruckerPragerXXVE_New(  DRUCKERPRAGERXXVE_PASSARGS  );
}

void _DruckerPragerXXVE_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
      DruckerPragerXXVE*            self           = (DruckerPragerXXVE*)rheology;
      FeVariable*           strainRateField;
      SwarmVariable*        swarmStrainRate;
      FeVariable*           pressureField;
      SwarmVariable*        swarmPressure;
      MaterialPointsSwarm*  materialPointsSwarm;
      IntegrationPointsSwarm* integrationSwarm;
      FiniteElementContext* context;

      /* required for initialization of plastic strain for strain weakening */

      double                  softeningStrain;
      double                  initialDamageFraction;
      double                  initialDamageWavenumber;
      double                  initialDamageWavenumberSinI;
      double                  initialDamageWavenumberCosI;
      double                  initialDamageWavenumberSinK;
      double                  initialDamageWavenumberCosK;
      double                  initialDamageFactor;
      long int                randomSeed;
      Stg_Shape*              initialStrainShape;

      /* Yielding and softening parameters */

      double cohesion, cohesionAfterSoftening, frictionCoefficient, frictionCoefficientAfterSoftening;
      double minimumViscosity, healingRate, minimumYieldStress;
      double alpha,alphaB,beta,kappa,kappaB;
      Bool   strainRateSoftening;

      /* Construct Parent */

      _Rheology_AssignFromXML( self, cf, data );
       Rheology_SetToNonLinear( rheology );

      strainRateField     = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"StrainRateField", FeVariable, False, data  );
      swarmStrainRate     = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"swarmStrainRate", SwarmVariable, False, data );
      pressureField       = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"PressureField", FeVariable, False, data  );
      swarmPressure       = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"swarmPressure", SwarmVariable, False, data );
      materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
      integrationSwarm    = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"IntegrationSwarm", IntegrationPointsSwarm, True, data  );
      context             = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  );
      
      Journal_Firewall( 
	    (strainRateField || self->swarmStrainRate ), 
	    Journal_Register( Error_Type, (Name)self->type  ), 
	    "\n Error in component type %s, name '%s'.\n Must specify a strainRateField OR a swarmStrainRate, but not both. \n", self->type, self->name ); 
      Journal_Firewall( 
	    (pressureField || self->swarmPressure ), 
	    Journal_Register( Error_Type, (Name)self->type  ), 
	    "\n Error in component type %s, name '%s'.\n Must specify a pressureField OR a swarmPressure, but not both. \n", self->type, self->name );
      
      /* required for initialization of plastic strain for strain weakening */
      softeningStrain             = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"softeningStrain", 1.0  );
      initialDamageFraction       = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageFraction", 0.0  );
      initialDamageWavenumber     = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageWavenumber", -1.0  );
      initialDamageWavenumberSinI = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageWavenumberSinI", -1.0  );
      initialDamageWavenumberCosI = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageWavenumberCosI", -1.0  );
      initialDamageWavenumberSinK = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageWavenumberSinK", -1.0  );
      initialDamageWavenumberCosK = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageWavenumberCosK", -1.0  );
      initialDamageFactor         = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialDamageFactor", 1.0 );
      randomSeed                  = (long int ) Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"randomSeed", 0  );
      initialStrainShape          = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"initialStrainShape", Stg_Shape, False, data  );

      /* stuff for viscoelasticity */

      self->viscoelastic_corr = False;
      self->veForceTerm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ViscoelasticForceTerm", ViscoelasticForceTerm, False, data  );
      if ( self->veForceTerm ) {
	       self->veRheology = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ViscoelasticRheology", ViscoelasticRheology, True, data  );
	       self->viscoelastic_corr  = True; /* flag that we have viscoelasticity on */
      } 
      else {
	       PetscPrintf(MPI_COMM_WORLD, "DruckerPragerXXVE: No viscoelastic force term was specified.\n");
	       PetscPrintf(MPI_COMM_WORLD, "DruckerPragerXXVE: No elastic contribution to the stress will be considered without this component.\n");
	       PetscPrintf(MPI_COMM_WORLD, "DruckerPragerXXVE: You should add <param name=\"ViscoelasticForceTerm\"> VE FORCE TERM COMPONENT NAME</param> to the component definition\n"); 
      }
   
 
      minimumYieldStress  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minimumYieldStress", 0.01  );
      frictionCoefficient = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficient", 0.0  );
      frictionCoefficientAfterSoftening = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"frictionCoefficientAfterSoftening", 0.0 );

      cohesion  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"cohesion", 1.0 );
      cohesionAfterSoftening = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"cohesionAfterSoftening", 0.01 );
      strainRateSoftening = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"strainRateSoftening", False );

      minimumViscosity = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minimumViscosity", 0.01  );

      healingRate = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"healingRate", 0.0  );

      alpha  = alphaB = 1.0/softeningStrain;
      beta = healingRate;
      if( frictionCoefficient > 0.0 ) {
	       kappaB = frictionCoefficientAfterSoftening/frictionCoefficient;
      }
      if(cohesion > 0.0) {
	       kappa = cohesionAfterSoftening/cohesion;
      }

      _DruckerPragerXXVE_Init( 
	         self,
	         materialPointsSwarm,
	         integrationSwarm,
	         context,
	         strainRateField,
	         swarmStrainRate,
	         pressureField,
	         swarmPressure,
	         cohesion,
	         frictionCoefficient,
	         minimumYieldStress,
	         minimumViscosity,
	         alpha,
	         beta,
	         kappa,
	         alphaB,
	         kappaB,
	         softeningStrain,
	         initialDamageFraction,
	         initialDamageWavenumber, 
	         initialDamageWavenumberSinI, 
	         initialDamageWavenumberCosI, 
	         initialDamageWavenumberSinK, 
	         initialDamageWavenumberCosK, 
	         initialDamageFactor,
	         randomSeed,
	         initialStrainShape );
}

void _DruckerPragerXXVE_Destroy( void* rheology, void* data ){
      DruckerPragerXXVE* self = (DruckerPragerXXVE*) rheology;
   
      if( self->strainRateField ) Stg_Component_Destroy( self->strainRateField, data, False );
      if( self->swarmStrainRate ) Stg_Component_Destroy( self->swarmStrainRate, data, False );
      if( self->pressureField ) Stg_Component_Destroy( self->pressureField, data, False );
      if( self->swarmPressure ) Stg_Component_Destroy( self->swarmPressure, data, False );
      if( self->plasticStrain ) Stg_Component_Destroy( self->plasticStrain, data, False );

      _Rheology_Destroy( self, data );
}

void _DruckerPragerXXVE_Build( void* rheology, void* data ){
      DruckerPragerXXVE* self = (DruckerPragerXXVE*) rheology;
   
      /* build parent */
      _Rheology_Build( self, data );
   
      if( self->strainRateField ) Stg_Component_Build( self->strainRateField, data, False );
      if( self->swarmStrainRate ) Stg_Component_Build( self->swarmStrainRate, data, False );
      if( self->pressureField ) Stg_Component_Build( self->pressureField, data, False );
      if( self->swarmPressure ) Stg_Component_Build( self->swarmPressure, data, False );
      if( self->plasticStrain ) Stg_Component_Build( self->plasticStrain, data, False );

}

void _DruckerPragerXXVE_Initialise( void* rheology, void* data ){
    DruckerPragerXXVE* self = (DruckerPragerXXVE*) rheology;
    Cell_LocalIndex          nElements           = self->materialPointsSwarm->cellLocalCount;
    MaterialPointsSwarm*     materialPointsSwarm = self->materialPointsSwarm;
    IntegrationPointsSwarm*  integrationSwarm    = self->integrationSwarm;
    DruckerPragerXXVE_ParticleExt*   particleExt;
    int                      elementIndex,   particleIndex;
    double                   plasticStrain;
    double*                  coord;

    /* Initialise parent */
    _Rheology_Initialise( self, data );
   
    if( self->strainRateField ) Stg_Component_Initialise( self->strainRateField, data, False );
    if( self->swarmStrainRate ) Stg_Component_Initialise( self->swarmStrainRate, data, False );
    if( self->pressureField ) Stg_Component_Initialise( self->pressureField, data, False );
    if( self->swarmPressure ) Stg_Component_Initialise( self->swarmPressure, data, False );
    Stg_Component_Initialise( self->materialPointsSwarm, data, False );
    Stg_Component_Initialise( self->integrationSwarm,    data, False );

    if ( self->context->loadSwarmsFromCheckpoint == False ) {
	Stg_Component_Initialise( self->plasticStrain, data, False );
	/* Initialise random number generator */
	srand( self->randomSeed );
	/* initialise the pstrain on the particles : for each particle in each element */
	for( elementIndex = 0; elementIndex < nElements; elementIndex++ ) {
	    int cellIndex  = CellLayout_MapElementIdToCellId( materialPointsSwarm->cellLayout, elementIndex );
	    int nParticles = materialPointsSwarm->cellParticleCountTbl[cellIndex];

	    for( particleIndex = 0; particleIndex < nParticles; particleIndex++ ) {
          IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleInCellAt( integrationSwarm, cellIndex, particleIndex );
          MaterialPoint*    matParticle = OneToOneMapper_GetMaterialPoint( integrationSwarm->mapper, integrationPoint, &(materialPointsSwarm) );
          Variable_SetValueFloat( self->plasticStrain->variable, particleIndex, 0.0 );
		
          /* Get Parameters From Material Swarm Extension */
      		particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, matParticle, self->particleExtHandle );
	
        	plasticStrain = 0.0;
		if ( rand() < RAND_MAX*self->initialDamageFraction ) {
		    coord = matParticle->coord;
		    if ( self->initialStrainShape && !Stg_Shape_IsCoordInside( self->initialStrainShape, coord ) ) {
	       		particleExt->pstrain = plasticStrain;
	       		continue;		   
        }

		    plasticStrain = self->initialDamageFactor * rand() * self->softeningStrain/RAND_MAX;
		    /* Modulate the initial weakening by a harmonic-squared function with wavenumber(s) specified by
		       the user. (Use old definition if new one is not set) */
		    if ( self->initialDamageWavenumber > 0.0 && self->initialDamageWavenumberCosI == -1.0 ) {				
			     plasticStrain *= pow(0.5+0.5*cos(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumber),2.0); 
        }
		    /* Alternate phase is appropriate for different bc's and choice of origin */	
		    if ( self->initialDamageWavenumberCosI > 0.0 ) {				
			     plasticStrain *= pow(0.5+0.5*cos(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumberCosI),2.0); 
        }
		    if ( self->initialDamageWavenumberSinI > 0.0 ) {				
			     plasticStrain *= pow(0.5+0.5*sin(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumberSinI),2.0); 
        }
		    if ( self->initialDamageWavenumberCosK > 0.0 ) {				
			     plasticStrain *= pow(0.5+0.5*cos(M_PI * coord[ K_AXIS ] * self->initialDamageWavenumberCosK),2.0); 
        }
		    if ( self->initialDamageWavenumberSinK > 0.0 ) {				
			     plasticStrain *= pow(0.5+0.5*sin(M_PI * coord[ K_AXIS ] * self->initialDamageWavenumberSinK),2.0); 
        }
		}

		particleExt->pstrain = plasticStrain;
		Variable_SetValueFloat( self->plasticStrain->variable, particleIndex, plasticStrain );
	    }/* each particle in element */
	   }/* each element */
    }/* if checkpointing false */   
  }

/* gets called by Rheology_Register_RunRheologies via the set function pointer (set via _DruckerPragerXXVE_DefaultNew) */

void _DruckerPragerXXVE_ModifyConstitutiveMatrix( 
      void*                                              yieldRheology, 
      ConstitutiveMatrix*                                constitutiveMatrix,
      MaterialPointsSwarm*                               materialPointsSwarm,
      Element_LocalIndex                                 lElement_I,
      MaterialPoint*                                     materialPoint,
      Coord                                              xi )   {

      DruckerPragerXXVE*  self          = (DruckerPragerXXVE*) yieldRheology;
      DruckerPragerXXVE_ParticleExt*  particleExt;
      Stokes_SLE*         SLE = (Stokes_SLE*)(constitutiveMatrix->sle);

      double      yieldCriterion;
      double      yieldIndicator; 
      double      c0,c1, gamma, dt, b0,b1;
      double      eta;
      double      strainInv;
      double      pressure; 
      SymmetricTensor strainRate; /* array of double with MAX_SYMMETRIC_TENSOR_COMPONENTS = 6 elements (StgDomain/Geometry/src/units.h) */
      double minYieldCriterion = self->minYieldCriterion;
      double minViscosity      = self->minViscosity;
      double alpha=self->alpha;
      double beta=self->beta;
      double kappa=self->kappa;
      double alphaB=self->alphaB;
      double kappaB=self->kappaB;

      /* Don't do anything if loading from checkpoint data and it's the first solve */
      if ( self->context->loadSwarmsFromCheckpoint == True && !constitutiveMatrix->previousSolutionExists ){
	       return;  
      }
      constitutiveMatrix->isDiagonal = True; /* just in case */

      /* the constitutiveMatrix data should be set 
         by the previous rheology in the rheology list in the system 
         which may be the ViscoelasticRheology */

      self->eta0 = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );

      /* Get Strain Rate */  /*  strainRate is filled with actual tensor values  */
      /* I haven't allowed for the existence of swarmStrainRate or swarmPressure in the DruckerPragerXXVE_IntegrateStrain function yet */
      /* because I don't have the constitutiveMatrix there */

      if( self->strainRateField ) {
	       FeVariable_InterpolateWithinElement( self->strainRateField, lElement_I, xi, strainRate );            
      }
      else {
	       SwarmVariable_ValueAt( self->swarmStrainRate, constitutiveMatrix->currentParticleIndex, strainRate ); 
      }
      if( self->pressureField ) {
	       FeVariable_InterpolateWithinElement( self->pressureField, lElement_I, xi, &pressure );             
      }
      else {
	     SwarmVariable_ValueAt( self->swarmPressure, constitutiveMatrix->currentParticleIndex, &pressure ); 
      }

      /*****************************************************************************************************************/

      /* apply viscoelastic correction */
      if ( self->viscoelastic_corr ){
        /* e -> e' = e + 1/(2*mu*dte)*tau_0; tau_0 = previous stress */

        /* get viscoelastic particle extension */
	       Viscoelastic_Particle* veExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->veForceTerm->particleExtHandle );
	       
         /* get viscoelastic rheology */
         ViscoelasticRheology* veRheology = self->veRheology;
	  
      	  double mu = veRheology->mu;
      	  double dt_e = veRheology->elasticTimeStep;
      	  double factor = 0.5 / (mu*dt_e);
      	  int dim=constitutiveMatrix->dim;
      	  int i;

	       for(i=0;i< dim*(dim+1)/2; i++){
	         strainRate[i] = strainRate[i] + factor * veExt->prevStress[i]; 
         }
      }

      /*****************************************************************************************************************/
      
      strainInv = SymmetricTensor_2ndInvariant( strainRate, constitutiveMatrix->dim );
      if(strainInv < 10*self->eps){/* skip if I2 ~ 0 */
	       return;
      }
      particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
      gamma = particleExt->pstrain;
      dt = constitutiveMatrix->context->dt;
      double eta0, dg;
      double Cterm, Bterm;
      /* the constitutiveMatrix data should be set by the previous rheology in the rheology list in the system */
      /* which may be the ViscoelasticRheology */
      eta0 = self->eta0;
      
      c0 = self->cohesion;
      c1 = kappa*c0;
      b0 = self->frictionCoeff;
      b1 = kappaB*b0;

      Cterm = (c0-c1)*exp(-gamma*alpha);
      Bterm = (b0-b1)*exp(-gamma*alphaB);
      
      yieldCriterion = Cterm + c1  + pressure*( Bterm + b1) ;
      eta = yieldCriterion/(2.0*strainInv);
      dg = ( 1.0 - (1.0+beta)*eta/eta0 )*strainInv*dt;/* Euler integrated delta term for current time-step y_(n+1) = y_n + F_n*dt */
      gamma += dg;
      if(gamma < 0.0){ gamma=0.0; }
      Cterm = (c0-c1)*exp(-gamma*alpha);
      Bterm = (b0-b1)*exp(-gamma*alphaB);

      yieldCriterion = Cterm + c1  + pressure*( Bterm + b1) ;/* up to date value */
      yieldIndicator = 2*eta0*strainInv;
      
      /*********************************************************************/
      int flag=0;
      if(yieldCriterion < minYieldCriterion) {
	  yieldCriterion = minYieldCriterion;
	  //pressure = (yieldCriterion - (Cterm + c1))/(Bterm + b1);
	  flag=1;
      }
      /*********************************************************************/

      /* Test to see if it has yielded */
      int dim = constitutiveMatrix->dim;
      int i,j;

      for(i=0;i< dim*(dim+1)/2; i++){
	  constitutiveMatrix->derivs[i] = 0.0;
      }

      if ( yieldIndicator > yieldCriterion ) {
	  double Idt, I2, I2_32;
	  int flagVisc=0;

	  Idt = strainInv*dt;
	  I2 = strainInv*strainInv;/* strain-rate invariant I2*/
	  I2_32 = I2*strainInv;   /* I2^(3/2) */

	  eta = yieldCriterion/(2.0*strainInv);
	  if( eta < minViscosity){
	      eta=minViscosity;
	      flagVisc=1;
	  }
	  
	  ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, eta );
	  
	  /* if we are here we want to modify the constitutiveMatrix->matrixData if we are building Jacobian for Newton-Krylov iteration */
	  if( constitutiveMatrix->sle && constitutiveMatrix->sle->nlFormJacobian && !flagVisc) {
	      double**   C  = constitutiveMatrix->matrixData;
	      double fac; /* common factor */
	      double fudge=sqrt(SLE->fnorm);
	      
	      /**************************************************************************************
		           fac = d viscosity / d I_2
		           and viscosity = DPXX viscosity
                           
                           strainInv = sqrt(I2) here
	      ***************************************************************************************/
		
	      double BB,CC,XX,AA,GG;
	      if(!flag){
		  BB = -Bterm*pressure*alphaB;
		  CC = -Cterm*alpha;
		  XX = 1.0 + dt*( BB + CC )*(1.0+beta)/(2.0*eta0);
		
		  /* if XX == 0.0 ever , then we cannot solve for d eta_y / d I_2 : can this happen? */
		  /*
		    Bterm and Cterm are always > 0 but the pressure could make things interesting.
		    If the background viscosity is small (making XX-1.0 larger in magnitude) then the velocity solution
		    should give higher velocities and hence a smaller dt via the Courant factor.
		  */
		
		  if( (XX)*(XX) < 1e-15 ){ fprintf(stderr,"Something has gone horribly wrong in %s: XX=0.0!\n", __func__); exit(1); }
		  
		  AA = dt*(CC+BB)*(1.0-(1.0+beta )*eta/eta0)/(4.0*I2);
	      }
	      GG = -(yieldCriterion)/(4.0*I2_32);
		  
	      //fac = AA + GG; /* first term is 0 if we set c+b*p = constant */
	      /* double check this for when we are at min yield.. */
	      if(!flag){
		  fac = 2.0*(AA+GG)/(XX); /* 2.0 * d_etavm/d_I2 */
	      } else { 
		  fac = GG; 	
	      }
	      /* The Jacobian addition to the C matrix */
	      for(i=0;i< dim*(dim+1)/2; i++){
		  for(j=0;j< dim*(dim+1)/2; j++){
		      C[i][j] += fac*strainRate[i]*strainRate[j];
		  }
	      }
	      for(i=0;i< dim*(dim+1)/2; i++){
		  /* if(self->context->timeStep > 1){ */
		  /*     if(C[i][i] < eta/100000000.0){ C[i][i] = eta/100000000.0; }/\* put a cap on how low the viscosity terms can go *\/ */
		  /* } */
		  /* else{ */
		  /*     if(C[i][i] < eta/100000.0){ C[i][i] = eta/100000.0; } */
		  /* } */
		  if(fudge < 1.0){
		      if(C[i][i] < eta*fudge*1e-3){ C[i][i] = eta*fudge*1e-3; }
		  }else{
		      if(C[i][i] < eta*0.1*fudge){ C[i][i] = eta*0.1*fudge; } 
		  }
		  
	      }
	      constitutiveMatrix->isDiagonal = False;
	      /* Now build the vector we are going to pass along to NonLinearGradientTermX */
	      GG = ( Bterm + b1 )/(2.0*strainInv);
	      fac = GG;  /* term is 0 if we set c+b*p = constant */
	      if(flag){ fac = 0.0; }
	      fac = 2.0*fac/XX;
	      for(i=0;i< dim*(dim+1)/2; i++){
		  constitutiveMatrix->derivs[i] = fac*strainRate[i];
	      }
		
	  }/*if constitutiveMatrix->sle && constitutiveMatrix->sle->nlFormJacobian */
      }/* if stress > criterion */
}
