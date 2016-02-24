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
#include "NonNewtonianAbs.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type NonNewtonianAbs_Type = "NonNewtonianAbs";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
NonNewtonianAbs* _NonNewtonianAbs_New(  NonNewtonianAbs_DEFARGS  ) 
{
	NonNewtonianAbs*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(NonNewtonianAbs) );
	self = (NonNewtonianAbs*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _NonNewtonianAbs_Init( NonNewtonianAbs* self, FeVariable* strainRateInvField, double initialViscosity, double strainRate, double stressExponent ) {

	self->strainRateInvField = strainRateInvField;
	self->initialViscosity = initialViscosity;
    self->defaultStrainRateInvariant = strainRate;
    self->stressExponent = stressExponent;

	Rheology_SetToNonLinear( self );
}

void _NonNewtonianAbs_Build( void* _self, void* data ){
	NonNewtonianAbs*  self = (NonNewtonianAbs*)_self;

	_Rheology_Build( self, data );
	
   Stg_Component_Build( self->strainRateInvField, data, False );
}

void _NonNewtonianAbs_Initialise( void* _self, void* data ){
	NonNewtonianAbs*  self = (NonNewtonianAbs*)_self;

   _Rheology_Initialise( self, data );

   Stg_Component_Initialise( self->strainRateInvField, data, False );
}

void _NonNewtonianAbs_Destroy( void* _self, void* data ){
	NonNewtonianAbs*  self = (NonNewtonianAbs*)_self;

   Stg_Component_Destroy( self->strainRateInvField, data, False );

	_Rheology_Destroy( self, data );
}


void* _NonNewtonianAbs_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(NonNewtonianAbs);
	Type                                                             type = NonNewtonianAbs_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _NonNewtonianAbs_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _NonNewtonianAbs_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _NonNewtonianAbs_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _NonNewtonianAbs_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _NonNewtonianAbs_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _NonNewtonianAbs_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _NonNewtonianAbs_New(  NonNewtonianAbs_PASSARGS  );
}

void _NonNewtonianAbs_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	NonNewtonianAbs*  self = (NonNewtonianAbs*)rheology;
	FeVariable*    strainRateInvField;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
	
	/* TODO: 'Keyfallback' soon to be deprecated/updated */
	strainRateInvField = Stg_ComponentFactory_ConstructByNameWithKeyFallback( cf, self->name, (Name)"StrainRateInvariantField", (Dictionary_Entry_Key)"StrainRateInvariantField", FeVariable, True, data  );
	/*strainRateInvField = Stg_ComponentFactory_ConstructByKey( cf, self->name,
				"StrainRateInvariantField", FeVariable, True);*/
				
				
				
	_NonNewtonianAbs_Init( 
			self,
			strainRateInvField,
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"initialViscosity", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"defaultStrainRateInvariant", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"stressExponent", 1.0 ) );
}

void _NonNewtonianAbs_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	NonNewtonianAbs*	              self              = (NonNewtonianAbs*) rheology;
	double                        strainRateInv;
	double                        viscosity;
	double                        n;

	/* On the first ever solve we use the default strain rate */
	if ( !constitutiveMatrix->previousSolutionExists ) {
	  strainRateInv = self->defaultStrainRateInvariant;
   } else {
     FeVariable_InterpolateWithinElement( self->strainRateInvField, lElement_I, xi, &strainRateInv );
   }

	if ( strainRateInv == 0 )
	  return;

	n = self->stressExponent;

	if ( n == 0.0 )
		return;

	/* Calculate New Viscosity */
	viscosity = self->initialViscosity * pow(strainRateInv, 1.0/n - 1.0);
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}


