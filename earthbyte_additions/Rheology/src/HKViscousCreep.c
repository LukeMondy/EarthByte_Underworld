#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "HKViscousCreep.h"

#include <float.h>
#include <assert.h>


/* Textual name of this class */
const Type HKViscousCreep_Type = "HKViscousCreep";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
HKViscousCreep* _HKViscousCreep_New(  HKViscousCreep_DEFARGS  ) 
{
   HKViscousCreep* self;

   /* Call private constructor of parent */
   assert( _sizeOfSelf >= sizeof(HKViscousCreep) );
   self = (HKViscousCreep*) _Rheology_New(  RHEOLOGY_PASSARGS  );
   
   return self;
}

void _HKViscousCreep_Init( HKViscousCreep* self, int srf, int temp, int pres, double n, double dsr, double A, int d, int wf, int mf, double p, double r, double alpha, double E, double V ) {

   self->strainRateInvTag = srf;
   self->temperatureLabel = temp;
   self->pressureLabel = pres;
   self->meltFractionLabel = mf;
   self->stressExponent = n;
   self->defaultStrainRateInvariant = dsr;
   self->preExponentialFactor = A;
   self->grainSizeLabel = d;
   self->grainSizeExponent = p;
   self->waterFugacityLabel = wf;
   self->waterFugacityExponent = r;
   self->alpha = alpha;
   self->activationEnergy = E;
   self->activationVolume = V;

   /* set to nonLinear rheology if n > 1 */
   if( n > 1 )
      Rheology_SetToNonLinear( self );
}


void _HKViscousCreep_Build( void* _self, void* data ){
   HKViscousCreep*  self = (HKViscousCreep*)_self;

   _Rheology_Build( self, data );
   
}

void _HKViscousCreep_Initialise( void* _self, void* data ){
   HKViscousCreep*  self = (HKViscousCreep*)_self;

   _Rheology_Initialise( self, data );
}

void _HKViscousCreep_Destroy( void* _self, void* data ){
   HKViscousCreep*  self = (HKViscousCreep*)_self;

   _Rheology_Destroy( self, data );
}


void* _HKViscousCreep_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                     _sizeOfSelf = sizeof(HKViscousCreep);
   Type                                                             type = HKViscousCreep_Type;
   Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
   Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
   Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
   Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _HKViscousCreep_DefaultNew;
   Stg_Component_ConstructFunction*                           _construct = _HKViscousCreep_AssignFromXML;
   Stg_Component_BuildFunction*                                   _build = _HKViscousCreep_Build;
   Stg_Component_InitialiseFunction*                         _initialise = _HKViscousCreep_Initialise;
   Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
   Stg_Component_DestroyFunction*                               _destroy = _HKViscousCreep_Destroy;
   Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _HKViscousCreep_ModifyConstitutiveMatrix;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*) _HKViscousCreep_New(  HKViscousCreep_PASSARGS  );
}


void _HKViscousCreep_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
   HKViscousCreep*  self = (HKViscousCreep*)rheology;
   PpcManager *mgr=NULL;
   int sr, temp, press;

   /* Construct Parent */
   _Rheology_AssignFromXML( self, cf, data );

   mgr = self->mgr;

   sr = PpcManager_GetField( self->mgr, cf, (Stg_Component*)self, "StrainRateInvariantField", False );
   temp = PpcManager_GetField( self->mgr, cf, (Stg_Component*)self, "TemperatureField", False );
   press = PpcManager_GetField( self->mgr, cf, (Stg_Component*)self, "PressureField", False );

   _HKViscousCreep_Init(
         self,
         sr,
         temp,
         press,
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StressExponent", 1.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"DefaultStrainRateInvariant", 1.0e-13 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"PreExponentialFactor", 1.0 ),
         PpcManager_GetPpcFromDict( mgr, cf, self->name, "GrainSize", "" ),
         PpcManager_GetPpcFromDict( mgr, cf, self->name, "WaterFugacity", "" ),
         PpcManager_GetPpcFromDict( mgr, cf, self->name, "MeltFraction", "" ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"GrainSizeExponent", 0.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"WaterFugacityExponent", 0.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MeltFractionFactor", 1.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"ActivationEnergy", 0.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"ActivationVolume", 0.0 ) );
}


void _HKViscousCreep_ModifyConstitutiveMatrix( void* _self, ConstitutiveMatrix* constitutiveMatrix, MaterialPointsSwarm* swarm, Element_LocalIndex lElement_I, MaterialPoint* materialPoint, Coord xi ) {
  HKViscousCreep* self = (HKViscousCreep*) _self;
  int err;
  double eII;
  double viscosity;
  double n, A, d, p, fH2O, r, a, F, E, V, temp, pres, R;
   /* stupid way to recover the integration point from the materialPoint
    * TODO: FIXCONSTITUTIVE use integration point only */
  IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleAt( constitutiveMatrix->integrationSwarm,
                                                                             constitutiveMatrix->currentParticleIndex );
  
  /** get parameters **/
  n = self->stressExponent;
  A = self->preExponentialFactor;
  
  /* get strain rate */
   if ( !constitutiveMatrix->previousSolutionExists ) {
      /* first solve uses default strain rate */
      eII = self->defaultStrainRateInvariant;
   } else {
      err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->strainRateInvTag, &eII );
      if( err ) { eII = self->defaultStrainRateInvariant; }
   }  

  /* get water */
  err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->waterFugacityLabel, &fH2O );
  if( err ) fH2O = 0;
  
  /* get melt */
  err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->meltFractionLabel, &F );
  if( err ) F = 0;
  
  /* get grain size */
  err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->grainSizeLabel, &d );
  if( err ) d = 0;

  p = self->grainSizeExponent;
  r = self->waterFugacityExponent;
  a = self->alpha;
  E = self->activationEnergy;
  V = self->activationVolume;
  /* get gas constant */
  PpcManager_GetConstant( self->mgr, "R", &R ); 
  
  /** Calculate New Viscosity **/

  viscosity = 0.5 * pow( A, -1.0/n );
  /* strain rate dependency */
  if( fabs( n - 1.0 ) > 1.0e-5 )
    viscosity *= pow( eII, 1.0/n - 1.0 );
  /* grain size dependency */
  if( p>DBL_EPSILON && d>DBL_EPSILON )
    viscosity *= pow( d, p/n );
  /* water dependency */
  if( r>DBL_EPSILON && fH2O>DBL_EPSILON )
    viscosity *= pow( fH2O, -r/n );
  /* melt dependency */
  if( F>DBL_EPSILON )
    viscosity *= exp( -a*F/n );

  if(  self->temperatureLabel > -1 ) {
     /* temperature and pressure dependency - arrhenius */

     /* get temperature */
     err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->temperatureLabel, &temp );
     assert( !err );

     /* get pressure */
     err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->pressureLabel, &pres );
     if( err ) pres = 0;

     viscosity *= exp( (E + pres*V) / (n*R*temp) );
  }
  if (viscosity == 0) {
    printf("pres %g, eII %g, A %g, n %f, temp %f", pres, eII, A, n, temp);
  }
  
  Journal_Firewall( !isnan(viscosity), Journal_Register( Error_Type, (Name)HKViscousCreep_Type ), "Viscosity is nan\n" );
  Journal_Firewall( viscosity!=0, Journal_Register( Error_Type, (Name)HKViscousCreep_Type ), "Viscosity is zero\n" );



  ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}


