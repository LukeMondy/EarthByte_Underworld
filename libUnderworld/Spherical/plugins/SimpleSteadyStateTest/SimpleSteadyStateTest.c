/**
	This plugin was created to easily switch out cartesian based algorithms for spherical based algorithms

	The advantages of using a plugin whilst developing the spherical code is that i don't have to spend too much time considering object design as I have access to all data structures in this plugin.
	**/


#include <mpi.h>
#include <float.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Spherical/Spherical.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

typedef struct
{
   __Codelet
   Vec prevT;
   Vec deltaT;
   /* THINGS TO SAVE */
} SimpleSteadyStatePlugin;

const Type SimpleSteadyStatePlugin_Type = "Spherical_SimpleSteadyStateTest";

void SimpleSteadyStatePlugin_CheckReached( UnderworldContext* context ){
   // get the advection diffusion data structure big hack to do by name
   AdvectionDiffusionSLE *sle = (AdvectionDiffusionSLE*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"EnergyEqn" );
   SimpleSteadyStatePlugin *self = (SimpleSteadyStatePlugin*)LiveComponentRegister_Get( context->CF->LCRegister, SimpleSteadyStatePlugin_Type );
   double infNorm;

   assert(self);
   assert(sle);

   // if first time step we create vectors prevT and deltaT. We also zero prevT.
   if( context->timeStep==0 || (context->timeStepSinceJobRestart<=1 && context->loadFromCheckPoint) ) {
      VecDuplicate( sle->phiVector->vector, &(self->prevT) );
      VecZeroEntries( self->prevT );
      VecDuplicate( sle->phiVector->vector, &(self->deltaT) );
      StgFEM_FrequentOutput_PrintValue( context, HUGE_VAL );
      return ;
   }

   // dT = T_[n] - T_[n-1]
   VecWAXPY( self->deltaT, -1, self->prevT, sle->phiVector->vector);
   VecNorm( self->deltaT, NORM_INFINITY, &infNorm);

   // copy prevT = current T
   VecCopy( sle->phiVector->vector, self->prevT );

   // print the norm value
   StgFEM_FrequentOutput_PrintValue( context, infNorm );

   if( infNorm < 1e-7 ) {
      /* If the norm is low enough then ask the context to stop running 
         at the conclusion of the current timestep. */
      printf("Inf norm of the temperature variation is %g - steady state - ending model\n", infNorm ); 
      context->gracefulQuit = True;

      // clean petsc vectors
      VecDestroy( self->prevT );
      VecDestroy( self->deltaT );
   }

}

void _SimpleSteadyStatePlugin_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )
{
   //SimpleSteadyStatePlugin* self = (SimpleSteadyStatePlugin*)component;
   UnderworldContext* context = NULL;

   context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data );

   StgFEM_FrequentOutput_PrintString( context, "inf_norm_deltaT" );
   ContextEP_Append( context, "Context_FrequentOutput", SimpleSteadyStatePlugin_CheckReached );
}


void* _SimpleSteadyStatePlugin_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(SimpleSteadyStatePlugin);
   Type                                                      type = SimpleSteadyStatePlugin_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SimpleSteadyStatePlugin_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _SimpleSteadyStatePlugin_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Codelet_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

   /* default value NON_GLOBAL */
   AllocationType nameAllocationType = NON_GLOBAL;

   return _Codelet_New( CODELET_PASSARGS );
}


Index Spherical_SimpleSteadyStateTest_Register( PluginsManager* pluginsManager )
{
   return PluginsManager_Submit( pluginsManager, SimpleSteadyStatePlugin_Type, (Name)"0", _SimpleSteadyStatePlugin_DefaultNew  );
}
