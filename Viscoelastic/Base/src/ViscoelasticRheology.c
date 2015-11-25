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

#include <string.h>
#include <assert.h> 

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type ViscoelasticRheology_Type = "ViscoelasticRheology";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
ViscoelasticRheology* _ViscoelasticRheology_New(  VISCOELASTICRHEOLOGY_DEFARGS  ) 
{
	ViscoelasticRheology*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(ViscoelasticRheology) );
	self = (ViscoelasticRheology*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _ViscoelasticRheology_Init(
		ViscoelasticRheology*                              self,
		FiniteElementContext*                              context,
		double                                             elasticTimeStep,
		double                                             mu )
{
	/* Assign Pointers */
	self->elasticTimeStep         = elasticTimeStep;
	self->mu                      = mu;
	
	EP_AppendClassHook( context->calcDtEP, _ViscoelasticRheology_Timestep, self );
}

void* _ViscoelasticRheology_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(ViscoelasticRheology);
	Type                                                             type = ViscoelasticRheology_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _ViscoelasticRheology_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _ViscoelasticRheology_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _ViscoelasticRheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _Rheology_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _ViscoelasticRheology_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _ViscoelasticRheology_New(  VISCOELASTICRHEOLOGY_PASSARGS  );
}

void _ViscoelasticRheology_AssignFromXML( void* viscoelasticity, Stg_ComponentFactory* cf, void* data ){
	ViscoelasticRheology*            self = (ViscoelasticRheology*)viscoelasticity;
	
	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );

	self->veForceTerm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ViscoelasticForceTerm", ViscoelasticForceTerm, False, data  );

	
	_ViscoelasticRheology_Init(
			self ,
			Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  ), 
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"elasticTimeStep", 1.0  ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"mu", 0.0 )  );
}

void _ViscoelasticRheology_Initialise( void* rheology, void* data ) {
	ViscoelasticRheology*                self                  = (ViscoelasticRheology*) rheology;
	
	_Rheology_Initialise( self, data );
}

void _ViscoelasticRheology_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi ) 
{
	ViscoelasticRheology*				self			= (ViscoelasticRheology*) rheology;
	Viscoelastic_Particle* 				veExt;

	double 								dt_e			= self->elasticTimeStep;
	double 								mu				= self->mu;
	double 								eta_eff;
	double 								eta;
	
	eta = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );

	/* Store the orignal viscosity for later use */
	veExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->veForceTerm->particleExtHandle );
	veExt->ParticleOriginalViscosity = eta;

	/* Load eta_eff into the constitutiveMatrix */
	eta_eff = eta * mu * dt_e / (eta + mu*dt_e);
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, eta_eff );
}

double _ViscoelasticRheology_Timestep( ViscoelasticRheology* self, void* context ) {
	/*there's no need to add a test here because the smallest timestep is chosen
	* somewhere else*/
	return self->elasticTimeStep/3.0;
}



