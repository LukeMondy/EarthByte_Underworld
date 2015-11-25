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


#ifndef __ExperimentalUnderworld_ViscoelasticForceTerm_h__
#define __ExperimentalUnderworld_ViscoelasticForceTerm_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type ViscoelasticForceTerm_Type;
	
	typedef struct {
		SymmetricTensor	ParticleStrainRate;
		SymmetricTensor	ParticleViscousStrainRate;
		double			ParticleViscousStrainRateInv;
		double			ParticleOriginalViscosity;
		SymmetricTensor	elasticStress;
		SymmetricTensor	viscousStress;
		SymmetricTensor	totalStress;
		SymmetricTensor stressRate;
		SymmetricTensor	prevStress; /* previous deviatoric bit */	
	    //double          prevPressure; /* for compressible fluid */
	} Viscoelastic_Particle;	

	/** class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __ViscoelasticForceTerm \
		/* Parent info */ \
		__ForceTerm \
		/* Virtual functions go here */ \
		/* Other Info */ \
		Dimension_Index			dim; \
		/** Necessary since the Integration/Material Swarm breakup */ \
		MaterialPointsSwarm*	materialPointsSwarm; \
		ConstitutiveMatrix*     constitutiveMatrix ; \
		ExtensionInfo_Index		particleExtHandle; \
		FeVariable*				strainRateField; \
		FeVariable*             pressureField; \
		Materials_Register*		materials_Register; \
		SwarmVariable*			ForceStress; \
		SwarmVariable*			Stress; \
		SwarmVariable*			ParticleStrainRate; \
		SwarmVariable*			ParticleViscousStrainRate; \
		SwarmVariable*			ParticleViscousStrainRateInv; \
		SwarmVariable*			ParticleOriginalViscosity; \
		Bool					largeDef; \
		Bool					stressSmoothing; \
		Bool					rk4_rotations; \
		Bool					trackingPtclFlagged; \
		JaumannRotator*			jaumannRotator; 

	struct ViscoelasticForceTerm { __ViscoelasticForceTerm };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define VISCOELASTICFORCETERM_DEFARGS \
                FORCETERM_DEFARGS

	#define VISCOELASTICFORCETERM_PASSARGS \
                FORCETERM_PASSARGS

	ViscoelasticForceTerm* _ViscoelasticForceTerm_New(  VISCOELASTICFORCETERM_DEFARGS  );
	
	/* 'Stg_Component' implementations */

	void* _ViscoelasticForceTerm_DefaultNew( Name name );
	void _ViscoelasticForceTerm_Init(                        
 			void*						forceTerm,
 			FiniteElementContext*		context,
 			ConstitutiveMatrix*       constitutiveMatrix,
 			FeVariable*				strainRateField,                                             
 			MaterialPointsSwarm*		materialPointsSwarm,                                         
 			JaumannRotator*			jaumannRotator,                                              
 			Materials_Register*		materials_Register );                                       
	void _ViscoelasticForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data );
	void _ViscoelasticForceTerm_Initialise( void* forceTerm, void* data );
	void _ViscoelasticForceTerm_Destroy( void* _self, void* data );
	void _ViscoelasticForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVector );

	/* Private Functions */
	void _ViscoelasticForceTerm_UpdateStress( void* forceTerm, void* data );

#endif

