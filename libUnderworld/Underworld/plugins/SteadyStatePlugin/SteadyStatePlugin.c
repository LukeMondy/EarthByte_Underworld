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
#include <assert.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Underworld/Vrms/Vrms.h>

#include "SteadyStatePlugin.h"

const Type Underworld_SteadyStatePlugin_Type = "Underworld_SteadyStatePlugin";

void Underworld_SteadyStatePlugin_CheckReached( UnderworldContext* context ){
	static int sscounter = 0;
	static double	VrmsPrev;
	double	VrmsNow;
	Dictionary*             dictionary         = context->dictionary;
	int       res          = Dictionary_GetInt( dictionary, (Dictionary_Entry_Key)"elementResI"  );
	char* filename;
	Stream* steadyStateVrmsInfo        = Journal_Register( Info_Type, (Name)"SteadyState_DataStream" );
	static int streamHasBeenSetup      = 0;
    Underworld_SteadyState* self;
	Underworld_Vrms* uwvrms;

	self = (Underworld_SteadyState*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Underworld_SteadyStatePlugin_Type );
    /* NB: it will be difficult to refactor this to handle other things than VRMS as the 
     * steady state variable of interest until we have a more integrated approach to frequent
     * output values. */
	uwvrms = (Underworld_Vrms*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Underworld_Vrms_Type );
	
	/* Output Stream Setup */
	if( !streamHasBeenSetup  ) {
		Stg_asprintf( &filename, "SteadyStateVrmsInfo_%d.dat", res  );
		Stream_RedirectFile_WithPrependedPath( steadyStateVrmsInfo, context->outputPath, filename );
		Stream_SetAutoFlush( steadyStateVrmsInfo, True );
		Memory_Free( filename );
		streamHasBeenSetup = 1;
	}
	
	VrmsNow = uwvrms->vrms;
	
	/*printf("vrms*tol is %g VrmsPrev - VrmsNow is %g sscounter is %d\n", VrmsPrev*tolerance, VrmsPrev - VrmsNow, sscounter);*/
	
	if(fabs(VrmsPrev - VrmsNow) < VrmsPrev*self->tolerance){
		sscounter = sscounter + 1;
		/*printf("vrms is %lf sscounter is %d counter is %d\n",VrmsNow,sscounter,counter);*/
	}
	else{
		sscounter = 0;
	}
	if(sscounter >= self->steadySteps){
		printf("we have reached steady state (%d steps with less than %g change in"
            " chosen observable.\n", self->steadySteps, self->tolerance);
		Journal_Printf( steadyStateVrmsInfo,"%d %8g\n",context->timeStep, VrmsNow );
		/* This will ask the context to stop running at the conclusion of
           the current timestep. */
        context->gracefulQuit = True;
	}
    else {
	    VrmsPrev = VrmsNow;
    }
}

void _Underworld_SteadyStatePlugin_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
    Underworld_SteadyState*   self = (Underworld_SteadyState*)component;
	AbstractContext*          context;

	context = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey( cf,
        self, (Dictionary_Entry_Key)"Context", UnderworldContext, True, data );
	self->tolerance = Dictionary_GetDouble_WithDefault( cf->rootDict,
        (Dictionary_Entry_Key)"tolerance", 0.0001  );
    self->steadySteps =  Dictionary_GetUnsignedInt_WithDefault( cf->rootDict,
        (Dictionary_Entry_Key)"steadySteps", 100 );    
	
	ContextEP_Append( context, "Context_FrequentOutput",
        Underworld_SteadyStatePlugin_CheckReached );
}


void* _Underworld_SteadyStatePlugin_DefaultNew( Name name ) {
    	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_SteadyState);
	Type                                                      type = Underworld_SteadyStatePlugin_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_SteadyStatePlugin_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_SteadyStatePlugin_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Initialise;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS   );
}


Index Underworld_SteadyStatePlugin_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_SteadyStatePlugin_Type, (Name)"0", _Underworld_SteadyStatePlugin_DefaultNew  );
}


