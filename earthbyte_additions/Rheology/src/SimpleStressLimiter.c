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
#include "SimpleStressLimiter.h"

#include <assert.h>
#include <float.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type SimpleStressLimiter_Type = "SimpleStressLimiter";

/* Public Constructor 
SimpleStressLimiter* SimpleStressLimiter_New(
      Name                  name,
      AbstractContext*      context,
      FeVariable*           temperatureField, 
      double                eta0,
      double                theta )
{
   SimpleStressLimiter* self = (SimpleStressLimiter*) _SimpleStressLimiter_DefaultNew( name );

   _Rheology_Init( self, (PICelleratorContext*)context );
   _SimpleStressLimiter_Init( self, temperatureField, eta0, theta );
   self->isConstructed = True;
   return self;
}*/

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
SimpleStressLimiter* _SimpleStressLimiter_New(  SimpleStressLimiter_DEFARGS  ) 
{
	SimpleStressLimiter*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(SimpleStressLimiter) );
	self = (SimpleStressLimiter*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _SimpleStressLimiter_Init( SimpleStressLimiter* self,  
                                  int                    srf, 
                                  double                 maxYieldStress, 
                                  double                 referenceStrainRate ) {
	self->strainRateInvTag = srf;
	self->maxYieldStress  = maxYieldStress;
    self->referenceStrainRate = referenceStrainRate;
}

void* _SimpleStressLimiter_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(SimpleStressLimiter);
	Type                                                             type = SimpleStressLimiter_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _SimpleStressLimiter_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _SimpleStressLimiter_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Rheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _SimpleStressLimiter_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _SimpleStressLimiter_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _SimpleStressLimiter_New(  SimpleStressLimiter_PASSARGS  );
}

void _SimpleStressLimiter_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	SimpleStressLimiter*  self = (SimpleStressLimiter*)rheology;
	int sr;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );

    Rheology_SetToNonLinear( self );
	
	sr = PpcManager_GetField( self->mgr, cf, (Stg_Component*)self, "StrainRateInvariantField", False );

	_SimpleStressLimiter_Init( 
			self, 
			sr,
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxYieldStress", 1.0  ),
            Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"referenceStrainRate", 1.0  )  );

    /* test if extension already exists */
   self->particleExtHandle = ExtensionManager_GetHandle( self->mgr->materialSwarm->particleExtensionMgr, (Name)"ClipStress" );
   if( self->particleExtHandle == (unsigned)-1 ) {
      MaterialPoint materialPoint;
      int *mechanismExtension;
      /* add the extension for a int */
      self->particleExtHandle = ExtensionManager_Add( self->mgr->materialSwarm->particleExtensionMgr, "ClipStress", sizeof(int) );
      mechanismExtension = (int*)ExtensionManager_Get( self->mgr->materialSwarm->particleExtensionMgr, &materialPoint, self->particleExtHandle );

      /* Setup extension variable to analyse dominate viscous mechanism: diffusion, dislocation, peierls */
      self->stressClip = Swarm_NewScalarVariable( self->mgr->materialSwarm, (Name)"ClippedStress", (ArithPointer) mechanismExtension - (ArithPointer) &materialPoint, Variable_DataType_Int );
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->stressClip->variable );
      LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*)self->stressClip );
   } else {
      self->stressClip = SwarmVariable_Register_GetByName( self->mgr->materialSwarm->swarmVariable_Register, "ClippedStress" );
   }
}

void _SimpleStressLimiter_Destroy( void* rheology, void* data ) {
	SimpleStressLimiter*          self               = (SimpleStressLimiter*) rheology;

	/* Destroy parent */
	_Rheology_Destroy( self, data );

}

void _SimpleStressLimiter_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
    SimpleStressLimiter*                 self              = (SimpleStressLimiter*) rheology;
    double      eII;
    double      viscosity;
    double      init_viscosity;
    double      current_stress;
    int         err;
    int *mechanismType=NULL;

    IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleAt( constitutiveMatrix->integrationSwarm, constitutiveMatrix->currentParticleIndex );

    viscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
    init_viscosity = viscosity;

    mechanismType = (int*)ExtensionManager_Get( swarm->particleExtensionMgr, materialPoint, self->particleExtHandle );
    assert(mechanismType);
    *mechanismType=-1;

    /* Calculate Parameters */  
    if ( !constitutiveMatrix->previousSolutionExists ) {
        /* first solve uses default strain rate */
        eII = self->referenceStrainRate;
    } else {
        err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->strainRateInvTag, &eII );
        if( err )
            eII = self->referenceStrainRate; 
    } 

    current_stress = 2.0 * init_viscosity * eII;

    if ( current_stress < self->maxYieldStress ) {
        *mechanismType=0;
        return;
    }
    else {
        *mechanismType=1;
        ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, self->maxYieldStress / (2*eII) );
    }
  
}


