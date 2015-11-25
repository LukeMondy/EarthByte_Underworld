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


#ifndef __ExperimentalUnderworld_JaumannRotator_h__
#define __ExperimentalUnderworld_JaumannRotator_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type JaumannRotator_Type;
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __JaumannRotator \
		/* Parent info */ \
 		__TimeIntegrand \
		/* Virtual functions go here */ \
		/* General Info */\
		/* Param passed in */ \
		FeVariable*                                         vorticityField;        \
		MaterialPointsSwarm*                                materialPointsSwarm;   \

	struct JaumannRotator { __JaumannRotator };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define JAUMANNROTATOR_DEFARGS \
                TIMEINTEGRAND_DEFARGS

	#define JAUMANNROTATOR_PASSARGS \
                TIMEINTEGRAND_PASSARGS

	JaumannRotator* _JaumannRotator_New(  JAUMANNROTATOR_DEFARGS  ) ;
	
	/* 'Stg_Component' implementations */
	void* _JaumannRotator_DefaultNew( Name name );
	void _JaumannRotator_AssignFromXML( void* JaumannRotator, Stg_ComponentFactory* cf, void* data );
	void _JaumannRotator_Build( void* JaumannRotator, void* data );
	void _JaumannRotator_Initialise( void* JaumannRotator, void* data );
	
	Bool _JaumannRotator_TimeDerivative( void* _JaumannRotator, Index lParticle_I, double* timeDeriv );
	void _JaumannRotator_Intermediate( void* JaumannRotator, Index lParticle_I );

	void JaumannRotator_UpdateVariables( void* timeIntegrator, JaumannRotator* self );

#endif

