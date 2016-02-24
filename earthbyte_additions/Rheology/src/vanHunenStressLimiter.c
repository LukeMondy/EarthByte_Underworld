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
#include "vanHunenStressLimiter.h"

#include <assert.h>
#include <float.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type vanHunenStressLimiter_Type = "vanHunenStressLimiter";

/* Public Constructor 
vanHunenStressLimiter* vanHunenStressLimiter_New(
      Name                  name,
      AbstractContext*      context,
      FeVariable*           temperatureField, 
      double                eta0,
      double                theta )
{
   vanHunenStressLimiter* self = (vanHunenStressLimiter*) _vanHunenStressLimiter_DefaultNew( name );

   _Rheology_Init( self, (PICelleratorContext*)context );
   _vanHunenStressLimiter_Init( self, temperatureField, eta0, theta );
   self->isConstructed = True;
   return self;
}*/

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
vanHunenStressLimiter* _vanHunenStressLimiter_New(  vanHunenStressLimiter_DEFARGS  ) 
{
	vanHunenStressLimiter*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(vanHunenStressLimiter) );
	self = (vanHunenStressLimiter*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _vanHunenStressLimiter_Init( vanHunenStressLimiter* self,  
                                  int                    srf, 
                                  double                 yieldStress, 
                                  double                 powerlawIndex, 
                                  double                 referenceStrainRate ) {
	self->strainRateInvTag = srf;
	self->yieldStress  = yieldStress;
	self->powerlawIndex = powerlawIndex;
  self->referenceStrainRate = referenceStrainRate;

  if ( powerlawIndex > 1 )
    Rheology_SetToNonLinear( self );

  Journal_Firewall( powerlawIndex != 0 , self->mgr->error_stream, "Error, in component %s.\nYou cannot have a powerlawIndex of 0.\n", self->name);

}

void* _vanHunenStressLimiter_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(vanHunenStressLimiter);
	Type                                                             type = vanHunenStressLimiter_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _vanHunenStressLimiter_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _vanHunenStressLimiter_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Rheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _vanHunenStressLimiter_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _vanHunenStressLimiter_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _vanHunenStressLimiter_New(  vanHunenStressLimiter_PASSARGS  );
}

void _vanHunenStressLimiter_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	vanHunenStressLimiter*  self = (vanHunenStressLimiter*)rheology;
	int sr;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
	
	sr = PpcManager_GetField( self->mgr, cf, (Stg_Component*)self, "StrainRateInvariantField", False );

	_vanHunenStressLimiter_Init( 
			self, 
			sr,
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"yieldStress", 1.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"powerlawIndex", 1.0 ),
      Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"referenceStrainRate", 1.0  )  );
}

void _vanHunenStressLimiter_Destroy( void* rheology, void* data ) {
	vanHunenStressLimiter*          self               = (vanHunenStressLimiter*) rheology;

	/* Destroy parent */
	_Rheology_Destroy( self, data );

}

void _vanHunenStressLimiter_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	vanHunenStressLimiter*                 self              = (vanHunenStressLimiter*) rheology;
	double      eII;
	double      viscosity;
  double      init_viscosity;
  int         err;
  IntegrationPoint* integrationPoint = (IntegrationPoint*) Swarm_ParticleAt( constitutiveMatrix->integrationSwarm,
                                                                             constitutiveMatrix->currentParticleIndex );
  
  viscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
  init_viscosity = viscosity;

	/* Calculate Parameters */  
  if ( !constitutiveMatrix->previousSolutionExists ) {
    /* first solve uses default strain rate */
    eII = self->referenceStrainRate;
  } else {
    err = PpcManager_Get( self->mgr, lElement_I, integrationPoint, self->strainRateInvTag, &eII );
    if( err ) { eII = self->referenceStrainRate; }
  } 

	/* Calculate New Viscosity */
	viscosity = self->yieldStress * pow(self->referenceStrainRate, (-1/self->powerlawIndex)) * pow(eII, (1/self->powerlawIndex)-1);
  if ( init_viscosity > viscosity ) {
	 ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
  }
  
}


