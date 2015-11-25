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
#include "ViscoelasticForceTerm.h"
#include "JaumannRotator.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type JaumannRotator_Type = "JaumannRotator";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
JaumannRotator* _JaumannRotator_New(  JAUMANNROTATOR_DEFARGS  ) 
{
	JaumannRotator*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree.
		At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(JaumannRotator) );
	self = (JaumannRotator*) _TimeIntegrand_New(  TIMEINTEGRAND_PASSARGS  );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _JaumannRotator_Init(
		JaumannRotator*                                    self,
		FeVariable*                                        vorticityField,
		MaterialPointsSwarm*                               materialPointsSwarm)
{
	/* Assign Values */
	self->vorticityField             = vorticityField;
	self->materialPointsSwarm        = materialPointsSwarm;
	/* NOTE: self->variable is set within _ViscoelasticForceTerm_Init routine */
	TimeIntegrator_AppendSetupEP( self->timeIntegrator,
		"JaumannRotator_UpdateVariables", JaumannRotator_UpdateVariables,  self->name, self );
		
}

void* _JaumannRotator_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(JaumannRotator);
	Type                                                       type = JaumannRotator_Type;
	Stg_Class_DeleteFunction*                               _delete = _TimeIntegrand_Delete;
	Stg_Class_PrintFunction*                                 _print = _TimeIntegrand_Print;
	Stg_Class_CopyFunction*                                   _copy = _TimeIntegrand_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = _JaumannRotator_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _JaumannRotator_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _JaumannRotator_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _JaumannRotator_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _TimeIntegrand_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _TimeIntegrand_Destroy;
	TimeIntegrand_CalculateTimeDerivFunction*  _calculateTimeDeriv = _JaumannRotator_TimeDerivative;
	TimeIntegrand_IntermediateFunction*              _intermediate = _JaumannRotator_Intermediate;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _JaumannRotator_New(  JAUMANNROTATOR_PASSARGS  );
}

void _JaumannRotator_AssignFromXML( void* _jaumannRotator, Stg_ComponentFactory* cf, void* data ){
	JaumannRotator* self = (JaumannRotator*) _jaumannRotator;
	MaterialPointsSwarm*    materialPointsSwarm;
	FeVariable*             vorticityField;
	
	/* Construct Parent */
	_TimeIntegrand_AssignFromXML( self, cf, data );
	
	vorticityField         = Stg_ComponentFactory_ConstructByNameWithKeyFallback( cf, self->name, (Name)"VorticityField", (Dictionary_Entry_Key)"VorticityField", FeVariable, True, data  );

	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
	
	_JaumannRotator_Init(
			self, 
			vorticityField,
			materialPointsSwarm);
			
			
}

void _JaumannRotator_Build( void* _jaumannRotator, void* data ) {
	JaumannRotator* self = (JaumannRotator*) _jaumannRotator;

	/* Build parent */
	_TimeIntegrand_Build( self, data );

}

void _JaumannRotator_Initialise( void* _jaumannRotator, void* data ) {
	JaumannRotator*   self = (JaumannRotator*) _jaumannRotator;
	
	/* Initialise Parent */
	_TimeIntegrand_Initialise( self, data );

	/* Update variables */
	Variable_Update( self->variable );

}

Bool _JaumannRotator_TimeDerivative( void* _jaumannRotator, Index lParticle_I, double* timeDeriv ) {
	JaumannRotator*          self = (JaumannRotator*) _jaumannRotator;
	MaterialPointsSwarm*     materialPointsSwarm = self->materialPointsSwarm;
	TensorArray              vorticity;
	Element_LocalIndex       lElement_I;
	MaterialPoint*           materialPoint = (MaterialPoint*) Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
	InterpolationResult      result;
	double*                  stress;
	double                   stressCalc[3][3], vorticityCalc[3][3] ;	
	double                   tauOmega_ji[3][3], tauOmega_ij[3][3] ;
	double                   k1[3][3] ;	
	int                      i, j, k;
	unsigned                 dim = materialPointsSwarm->dim;
	
					    
	stress = Variable_GetPtrDouble( self->variable, lParticle_I );

	lElement_I  = materialPoint->owningCell;
	
	result = FieldVariable_InterpolateValueAt( self->vorticityField, materialPoint->coord, vorticity );

	if ( result == OTHER_PROC || result == OUTSIDE_GLOBAL || isinf(vorticity[0]) || isinf(vorticity[1]) || 
		( dim == 3 && isinf(vorticity[2]) ) ) 
	{
		
		Journal_Printf( Journal_Register( Error_Type, (Name)self->type  ),
			"Error in func '%s' for particle with index %u.\n\tPosition (%g, %g, %g)\n\tvorticity here is (%g, %g, %g)."
			"\n\tInterpolation result is %s.\n",
			__func__, lParticle_I, materialPoint->coord[0], materialPoint->coord[1], materialPoint->coord[2], 
			vorticity[0], vorticity[1], ( dim == 3 ? vorticity[2] : 0.0 ),
			InterpolationResultToStringMap[result]  );
		return False;	
			
	}

	if ( dim == 2 ) {
		stressCalc[0][0] = stress[0];        stressCalc[0][1] = stress[2];
		stressCalc[1][0] = stress[2];        stressCalc[1][1] = stress[1];
		
		vorticityCalc[0][0] = vorticity[0];  vorticityCalc[0][1] = vorticity[1];
		vorticityCalc[1][0] = vorticity[2];  vorticityCalc[1][1] = vorticity[3];
	}
	else {
		stressCalc[0][0] = stress[0];       stressCalc[0][1] = stress[3];		stressCalc[0][2] = stress[4];
		stressCalc[1][0] = stress[3];       stressCalc[1][1] = stress[1];		stressCalc[1][2] = stress[5];
		stressCalc[2][0] = stress[4];       stressCalc[2][1] = stress[5];		stressCalc[2][2] = stress[2];
		
		vorticityCalc[0][0] = vorticity[0];		vorticityCalc[0][1] = vorticity[1];		vorticityCalc[0][2] = vorticity[2];
		vorticityCalc[1][0] = vorticity[3];		vorticityCalc[1][1] = vorticity[4];		vorticityCalc[1][2] = vorticity[5];
		vorticityCalc[2][0] = vorticity[6];		vorticityCalc[2][1] = vorticity[7];		vorticityCalc[2][2] = vorticity[8];
	}

	/* k1 - stress rate values from previous timestep */
	for ( i = 0 ; i < dim ; i++ ) {
		for ( j = 0 ; j < dim ; j++ ) {
			tauOmega_ji[i][j] = 0.0;
			tauOmega_ij[i][j] = 0.0;
			for ( k = 0 ; k < dim ; k++ ) {
				tauOmega_ij[i][j] = tauOmega_ij[i][j] + stressCalc[i][k] * vorticityCalc[k][j] ;
				tauOmega_ji[i][j] = tauOmega_ji[i][j] + stressCalc[j][k] * vorticityCalc[k][i] ;
			}
		}
	}
	for ( i = 0 ; i < dim ; i++ )	{
		for ( j = 0 ; j < dim ; j++ ) { 
			k1[i][j] = tauOmega_ij[i][j] + tauOmega_ji[i][j] ;
		}
	}

	if ( dim == 2 ) {
		timeDeriv[0] = -k1[0][0];
		timeDeriv[1] = -k1[1][1];
		timeDeriv[2] = -k1[0][1];
	}
	else {
		timeDeriv[0] = -k1[0][0];		timeDeriv[3] = -k1[0][1];		timeDeriv[4] = -k1[0][2];
		timeDeriv[1] = -k1[1][1];		timeDeriv[5] = -k1[1][2];		timeDeriv[2] = -k1[2][2];
	}	
	
	return True;
}

/* This function is called after each of the time integration steps */ 
void _JaumannRotator_Intermediate( void* jaumannRotator, Index lParticle_I ) {
}

/* Update these variables in case the swarm has been reallocated */
void JaumannRotator_UpdateVariables( void* timeIntegrator, JaumannRotator* self ) {
	Variable_Update( self->variable );
	Variable_Update( self->materialPointsSwarm->particleCoordVariable->variable );
	Variable_Update( self->materialPointsSwarm->owningCellVariable->variable );
}
	


