#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <string.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "Ppc_VonMisesFlowRule.h"


/* Textual name of this class */
const Type Ppc_VonMisesFlowRule_Type = "Ppc_VonMisesFlowRule";


void _Ppc_VonMisesFlowRule_Init( void* _self, 
      int          viscosityTag,
      int          strainRateInvTag,
      int          yieldCritTag ) {
   Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;
   MaterialPointsSwarm* ms = self->manager->materialSwarm;
   Index handle=-1;
   ArithPointer offset; 

   self->viscosityTag = viscosityTag;
   self->strainRateInvTag = strainRateInvTag;
   self->yieldCritTag = yieldCritTag;

   assert ( ms ); 
          
   /* See if the YieldRheology Type has already added an extension to the particle 
   * If handle is given a value of '-1' - that means that no extension has been added
   * with the YieldRheology_Type
   * We should then add the extension */
          
   handle = ExtensionManager_GetHandle( ms->particleExtensionMgr, "HasYieldedExtension" );
       
   if ( handle == (ExtensionInfo_Index) -1  ) {
      handle = ExtensionManager_Add( ms->particleExtensionMgr, "HasYieldedExtension", sizeof(Index));
       
      /* Adding variable for plotting purpose */
      offset = (ArithPointer) ExtensionManager_Get( ms->particleExtensionMgr, NULL, handle );

      self->hasYieldedVariable = Swarm_NewScalarVariable( ms, (Name)"HasYielded", offset, Variable_DataType_Char  ); 
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->hasYieldedVariable->variable ); 
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->hasYieldedVariable ); 
   }
  
  self->hasYieldedIndex = handle;
}
   
/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_VonMisesFlowRule* _Ppc_VonMisesFlowRule_New(  PPC_OPERATION_DEFARGS  )
{
	Ppc_VonMisesFlowRule* self;
	
	assert( _sizeOfSelf >= sizeof(Ppc_VonMisesFlowRule) );
	nameAllocationType = NON_GLOBAL;
	self = (Ppc_VonMisesFlowRule*) _Ppc_New(  PPC_PASSARGS  );	
	self->_get = _get;
	return self;
}


void* _Ppc_VonMisesFlowRule_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Ppc_VonMisesFlowRule);
	Type                                                      type = Ppc_VonMisesFlowRule_Type;
	Stg_Class_DeleteFunction*                              _delete = _Ppc_VonMisesFlowRule_Delete;
	Stg_Class_PrintFunction*                                _print = _Ppc_VonMisesFlowRule_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_VonMisesFlowRule_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Ppc_VonMisesFlowRule_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Ppc_VonMisesFlowRule_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Ppc_VonMisesFlowRule_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Ppc_VonMisesFlowRule_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Ppc_VonMisesFlowRule_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_VonMisesFlowRule_Get;

	return (void*)_Ppc_VonMisesFlowRule_New(  PPC_OPERATION_PASSARGS  );
}


void _Ppc_VonMisesFlowRule_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;
  
  /* Construct parent */
  _Ppc_AssignFromXML( self, cf, data );

   self->sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", SystemLinearEquations, True, data );

  _Ppc_VonMisesFlowRule_Init(
      self,
      PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"viscosity", "" ),
      PpcManager_GetField( self->manager, cf, self, (Dictionary_Entry_Key)"strainrateinv", True ),
      PpcManager_GetPpcFromDict( self->manager, cf, self->name, (Dictionary_Entry_Key)"yieldcriterion", "" ));
}


void _Ppc_VonMisesFlowRule_Build( void* _self, void* data ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;

  /* Build parent */
  _Ppc_Build( self, data );
  Stg_Component_Build( self->sle, data, False );

  // code to set to non linear solves always
  SystemLinearEquations_SetToNonLinear( self->sle );
}

void _Ppc_VonMisesFlowRule_Initialise( void* _self, void* data ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;
  SwarmVariable* hasYieldVar=NULL; // no need to store on class

  /* Initialize parent */
  _Ppc_Initialise( self, data );
  Stg_Component_Initialise( self->sle, data, False );

  /* probably should be a in another location because other yield mechanicsm should use it */
  if ( self->context->loadSwarmsFromCheckpoint == False ) { 
     MaterialPointsSwarm* ms=self->manager->materialSwarm;
     MaterialPoint *mp=NULL;
     int p_i, particleLocalCount, *ext; 
     particleLocalCount = ms->particleLocalCount;
     // note the original swarm variable would have been initialised previously my the ppcManager_initialise()
     for( p_i=0; p_i<particleLocalCount; p_i++ ) {
         mp = (MaterialPoint*)Swarm_ParticleAt( ms, p_i );
         ext = ExtensionManager_Get( ms->particleExtensionMgr, mp, self->hasYieldedIndex );
         *ext = 0;
     }
   }	
}

void _Ppc_VonMisesFlowRule_Delete( void* _self ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;

  /* Delete parent */
  _Ppc_Delete( self );
}

void _Ppc_VonMisesFlowRule_Print( void* _self, Stream* stream ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;

  /* Print parent */
  _Ppc_Print( self, stream );
}

void _Ppc_VonMisesFlowRule_Execute( void* _self, void* data ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self;

  /* Execute parent */
  _Ppc_Execute( self, data );
}

void _Ppc_VonMisesFlowRule_Destroy( void* _self, void* data ) {
  Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*)_self; 
  
  Stg_Component_Destroy( self->sle, data, False );

  /* Destroy parent */
  _Ppc_Destroy( self, data );
}


/* 
 * Public functions 
 *
 */

void SetHasYieldedExtension( Ppc_VonMisesFlowRule* self, IntegrationPoint* particle, int hasYielded ) {
  int* ext = (int*)_OneToOneMapper_GetExtensionOn(((IntegrationPointsSwarm*)self->manager->integrationSwarm)->mapper, particle, self->hasYieldedIndex );
  assert(ext);
  *ext = hasYielded;
}

int _Ppc_VonMisesFlowRule_Get( void* _self, Element_LocalIndex lElement_I, IntegrationPoint* particle, double* result ) {
   Ppc_VonMisesFlowRule* self = (Ppc_VonMisesFlowRule*) _self;
   PpcManager* mgr = self->manager;
   double eta, eII, tau_y, tauII;
   int hasYielded = 0;
   int err;

   err = PpcManager_Get( mgr, lElement_I, particle, self->viscosityTag, &eta);
   err = PpcManager_Get( mgr, lElement_I, particle, self->strainRateInvTag, &eII);
   err = PpcManager_Get( mgr, lElement_I, particle, self->yieldCritTag, &tau_y);

   tauII = 2*eta*eII; // should this be a ppc function, probably yes

   // ugly to check for loadFromCheckPoint here
   if( tauII >= tau_y && (self->sle->hasExecuted || self->context->loadFromCheckPoint) ) {
      // plastic time
      /* Is there strain weakening */ 
      /*
      if( self->strainRateWeakening ) { 
         err = PpcManager_Get( mgr, lElement_I, particle, self->strainRateWeakening, & ); 
         assert(!err);
         viscosity = 2.0 * eta * tau_y * tau_y / ( 
      } else {
      */
      eta = tau_y / (2*eII);
      hasYielded=1;
   }

   *result=eta;

   SetHasYieldedExtension( self, particle, hasYielded );
      
   return 0;
}
