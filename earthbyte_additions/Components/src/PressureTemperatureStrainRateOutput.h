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


#ifndef __Underworld_Utils_PressureTemperatureStrainRateOutput_h__
#define __Underworld_Utils_PressureTemperatureStrainRateOutput_h__

	/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type PressureTemperatureStrainRateOutput_Type;

	/* PressureTemperatureStrainRateOutput information */
	#define __PressureTemperatureStrainRateOutput \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__SwarmOutput \
		/* Virtual Info */\
		FeVariable*                                       pressureField;                \
		FeVariable*                                       temperatureField;             \
		FeVariable*                                       strainRateField;  

	struct PressureTemperatureStrainRateOutput { __PressureTemperatureStrainRateOutput };
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PressureTemperatureStrainRateOutput_DEFARGS \
                SWARMOUTPUT_DEFARGS

	#define PressureTemperatureStrainRateOutput_PASSARGS \
                SWARMOUTPUT_PASSARGS

	PressureTemperatureStrainRateOutput* _PressureTemperatureStrainRateOutput_New(  PressureTemperatureStrainRateOutput_DEFARGS  );

	/* Stg_Class_Delete PressureTemperatureStrainRateOutput implementation */
	void _PressureTemperatureStrainRateOutput_Delete( void* swarmOutput );
	void _PressureTemperatureStrainRateOutput_Print( void* swarmOutput, Stream* stream );
	#define PressureTemperatureStrainRateOutput_Copy( self ) \
		(PressureTemperatureStrainRateOutput*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define PressureTemperatureStrainRateOutput_DeepCopy( self ) \
		(PressureTemperatureStrainRateOutput*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _PressureTemperatureStrainRateOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _PressureTemperatureStrainRateOutput_DefaultNew( Name name ) ;
void _PressureTemperatureStrainRateOutput_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _PressureTemperatureStrainRateOutput_Build( void* swarmOutput, void* data ) ;
	void _PressureTemperatureStrainRateOutput_Initialise( void* swarmOutput, void* data ) ;
	void _PressureTemperatureStrainRateOutput_Execute( void* swarmOutput, void* data );
	void _PressureTemperatureStrainRateOutput_Destroy( void* swarmOutput, void* data ) ;
	
	void _PressureTemperatureStrainRateOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
	void _PressureTemperatureStrainRateOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context );
		
	/*---------------------------------------------------------------------------------------------------------------------
	** Private functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Entry Point Hooks
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/
	
#endif 

