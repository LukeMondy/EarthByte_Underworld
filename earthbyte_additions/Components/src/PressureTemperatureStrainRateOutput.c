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

#include "types.h"
#include "PressureTemperatureStrainRateOutput.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type PressureTemperatureStrainRateOutput_Type = "PressureTemperatureStrainRateOutput";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
PressureTemperatureStrainRateOutput* PressureTemperatureStrainRateOutput_New(
		Name                                  name,
		TimeIntegrator*                       timeIntegrator,
		FeVariable*                           velocityField )
{
	PressureTemperatureStrainRateOutput* self = (PressureTemperatureStrainRateOutput*) _PressureTemperatureStrainRateOutput_DefaultNew( name );

	/* 	PressureTemperatureStrainRateOutput_InitAll */
	abort();

	return self;
}

PressureTemperatureStrainRateOutput* _PressureTemperatureStrainRateOutput_New(  PressureTemperatureStrainRateOutput_DEFARGS  )
{
	PressureTemperatureStrainRateOutput* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(PressureTemperatureStrainRateOutput) );
	self = (PressureTemperatureStrainRateOutput*)_SwarmOutput_New(  SWARMOUTPUT_PASSARGS  );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _PressureTemperatureStrainRateOutput_Init( 
		void*                                 swarmOutput,
		FeVariable*                           pressureField,
		FeVariable*                           temperatureField,
		FeVariable*                           strainRateField )
{
	PressureTemperatureStrainRateOutput*   self                = (PressureTemperatureStrainRateOutput*)swarmOutput;

	self->pressureField      = pressureField;
	self->temperatureField   = temperatureField;
	self->strainRateField    = strainRateField;
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PressureTemperatureStrainRateOutput_Delete( void* swarmOutput ) {
	PressureTemperatureStrainRateOutput* self = (PressureTemperatureStrainRateOutput*)swarmOutput;

	/* Delete parent */
	_SwarmOutput_Delete( self );
}


void _PressureTemperatureStrainRateOutput_Print( void* swarmOutput, Stream* stream ) {
	PressureTemperatureStrainRateOutput* self = (PressureTemperatureStrainRateOutput*)swarmOutput;
	
	/* Print parent */
	_SwarmOutput_Print( self, stream );
}

void* _PressureTemperatureStrainRateOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*)swarmOutput;
	PressureTemperatureStrainRateOutput*	newPressureTemperatureStrainRateOutput;
	
	newPressureTemperatureStrainRateOutput = (PressureTemperatureStrainRateOutput*)_SwarmOutput_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newPressureTemperatureStrainRateOutput;
}

void* _PressureTemperatureStrainRateOutput_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(PressureTemperatureStrainRateOutput);
	Type                                                      type = PressureTemperatureStrainRateOutput_Type;
	Stg_Class_DeleteFunction*                              _delete = _PressureTemperatureStrainRateOutput_Delete;
	Stg_Class_PrintFunction*                                _print = _PressureTemperatureStrainRateOutput_Print;
	Stg_Class_CopyFunction*                                  _copy = _PressureTemperatureStrainRateOutput_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _PressureTemperatureStrainRateOutput_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _PressureTemperatureStrainRateOutput_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _PressureTemperatureStrainRateOutput_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _PressureTemperatureStrainRateOutput_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _PressureTemperatureStrainRateOutput_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _PressureTemperatureStrainRateOutput_Destroy;
	SwarmOutput_PrintHeaderFunction*                  _printHeader = _PressureTemperatureStrainRateOutput_PrintHeader;
	SwarmOutput_PrintDataFunction*                      _printData = _PressureTemperatureStrainRateOutput_PrintData;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _PressureTemperatureStrainRateOutput_New(  PressureTemperatureStrainRateOutput_PASSARGS  );
}


void _PressureTemperatureStrainRateOutput_AssignFromXML( void* swarmOutput, Stg_ComponentFactory* cf, void* data ) {
	PressureTemperatureStrainRateOutput*  self          = (PressureTemperatureStrainRateOutput*) swarmOutput;
	FeVariable*                 pressureField;
	FeVariable*                 temperatureField;
	FeVariable*                 strainRateField;

	_SwarmOutput_AssignFromXML( self, cf, data );

	pressureField    = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"PressureField", FeVariable, True, data  ) ;
	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"TemperatureField", FeVariable, True, data  ) ;
    strainRateField  = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"StrainRateField", FeVariable, True, data  ) ;
	_PressureTemperatureStrainRateOutput_Init( self, pressureField, temperatureField, strainRateField );

}

void _PressureTemperatureStrainRateOutput_Build( void* swarmOutput, void* data ) {
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*) swarmOutput;

	Stg_Component_Build( self->pressureField, data, False );
	Stg_Component_Build( self->temperatureField, data, False );
	Stg_Component_Build( self->strainRateField, data, False );

	_SwarmOutput_Build( self, data );
}
void _PressureTemperatureStrainRateOutput_Initialise( void* swarmOutput, void* data ) {
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*) swarmOutput;

	Stg_Component_Initialise( self->pressureField, data, False );
	Stg_Component_Initialise( self->temperatureField, data, False );
	Stg_Component_Initialise( self->strainRateField, data, False );
	
	_SwarmOutput_Initialise( self, data );
}
void _PressureTemperatureStrainRateOutput_Execute( void* swarmOutput, void* data ) {
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*)swarmOutput;
	
	_SwarmOutput_Execute( self, data );
}
void _PressureTemperatureStrainRateOutput_Destroy( void* swarmOutput, void* data ) {
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*)swarmOutput;

	Stg_Component_Destroy( self->pressureField, data, False );
	Stg_Component_Destroy( self->temperatureField, data, False );
	Stg_Component_Destroy( self->strainRateField, data, False );
	
	_SwarmOutput_Destroy( self, data );
}

void _PressureTemperatureStrainRateOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	PressureTemperatureStrainRateOutput*	self = (PressureTemperatureStrainRateOutput*)swarmOutput;
	
	_SwarmOutput_PrintHeader( self, stream, lParticle_I, context );
	
	SwarmOutput_PrintString( self, stream, "Pressure" );
	SwarmOutput_PrintString( self, stream, "Temperature" );
	SwarmOutput_PrintString( self, stream, "StrainRate" );
}

void _PressureTemperatureStrainRateOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	PressureTemperatureStrainRateOutput*	self                = (PressureTemperatureStrainRateOutput*)swarmOutput;
	Swarm*                      swarm               = self->swarm;
	GlobalParticle*           particle            = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
	double*                     coord               = particle->coord;
	double                      pressure;
	double                      temperature;
    double                      strainRate;

	Journal_Firewall(
		swarm->particleLayout->coordSystem == GlobalCoordSystem,
		Journal_MyStream( Error_Type, self ),
		"Swarm is not using global coord system! Modify this code to use both systems\n" );

	_SwarmOutput_PrintData( self, stream, lParticle_I, context );

	FieldVariable_InterpolateValueAt( self->pressureField,    coord, &pressure );
	FieldVariable_InterpolateValueAt( self->temperatureField, coord, &temperature );
	FieldVariable_InterpolateValueAt( self->strainRateField, coord, &strainRate );
	
	SwarmOutput_PrintValue( self, stream, pressure );
	SwarmOutput_PrintValue( self, stream, temperature );
	SwarmOutput_PrintValue( self, stream, strainRate );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
/*---------------------------------------------------------------------------------------------------------------------
** Entry Point Hooks
*/

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/



