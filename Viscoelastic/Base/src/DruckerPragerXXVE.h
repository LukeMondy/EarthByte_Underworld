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


#ifndef __Underworld_Rheology_DruckerPragerXXVE_h__
#define __Underworld_Rheology_DruckerPragerXXVE_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type DruckerPragerXXVE_Type;
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __DruckerPragerXXVE \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Rheology \
		/* Virtual functions go here */ \
		/* Material Parameters */\
		FeVariable*                                         strainRateField;                      \
		SwarmVariable*                                      swarmStrainRate;                      \
		FeVariable*                                         pressureField;                        \
		SwarmVariable*                                      swarmPressure;                        \
		double                                              cohesion;                             \
		double                                              frictionCoeff;                        \
        SymmetricTensor                                     stress;									\
        double                                              stressInv;								\
		ExtensionInfo_Index                                 particleExtHandle;						\
		MaterialPointsSwarm*                                materialPointsSwarm;\
		IntegrationPointsSwarm*                             integrationSwarm;\
		double												eta0;\
		double												minYieldCriterion;\
		double												minViscosity;\
		double												dt;\
		double												alpha;\
		double												beta;\
		double												kappa;\
		double												alphaB;\
		double												betaB;\
		double												kappaB;\
		double												eps;   \
		int   												epsI;  \
		double												softeningStrain;                   \
		double												initialDamageFraction;             \
		double												initialDamageWavenumber;           \
		double												initialDamageWavenumberSinI;           \
		double												initialDamageWavenumberCosI;           \
		double												initialDamageWavenumberSinK;           \
		double												initialDamageWavenumberCosK;           \
		double												initialDamageFactor;               \
		Stg_Shape*											initialStrainShape;                \
		long int  											randomSeed;  \
		SwarmVariable* 										plasticStrain; \
        SwarmVariable* 										tensileFailure; \
		/* Param passed in */																	\
		int                                                 pressureTag; \
		double                                              minimumYieldStress;	\
		double                                              frictionCoefficientAfterSoftening; \
		Bool                                                viscoelastic_corr; \
		ViscoelasticForceTerm* 								veForceTerm;		       \
		ViscoelasticRheology* 								veRheology;			\
		double 												curFrictionCoef;

	struct DruckerPragerXXVE { __DruckerPragerXXVE };
 
    typedef struct  {
	  	double pstrain;
	  	int pstrainIntegrated;
	} DruckerPragerXXVE_ParticleExt;
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define DRUCKERPRAGERXXVE_DEFARGS \
                RHEOLOGY_DEFARGS

	#define DRUCKERPRAGERXXVE_PASSARGS \
                RHEOLOGY_PASSARGS

	DruckerPragerXXVE* _DruckerPragerXXVE_New(  DRUCKERPRAGERXXVE_DEFARGS  );
	
	/* 'Stg_Component' implementations */
        
    void* _DruckerPragerXXVE_DefaultNew( Name name ) ;
	void  _DruckerPragerXXVE_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
    void  _DruckerPragerXXVE_Init( 
    			DruckerPragerXXVE* 			self, 
				MaterialPointsSwarm* 		materialPointsSwarm, 
				IntegrationPointsSwarm*		integrationSwarm, 
				FiniteElementContext*		context, 
				FeVariable*					strainRateField, 
				SwarmVariable* 				swarmStrainRate, 
				FeVariable*					pressureField, 
				SwarmVariable* 				swarmPressure, 
				double						cohesion, 
				double						frictionCoeff,
				double						minYieldCriterion,
				double						minViscosity,	   
				double						alpha, 
				double						beta, 
				double						kappa,
				double						alphaB, 
				double						kappaB,
				double						softeningStrain,
				double						initialDamageFraction,
				double						initialDamageWavenumber,
				double						initialDamageWavenumberSinI,
				double						initialDamageWavenumberCosI,
				double						initialDamageWavenumberSinK,
				double						initialDamageWavenumberCosK,
				double						initialDamageFactor,
				long int   					randomSeed,
				Stg_Shape* 					initialStrainShape
			);

        void _DruckerPragerXXVE_UpdateDrawingVariables( void* rheology );
        void  DruckerPragerXXVE_IntegrateStrain( void* rheology );
        void  DruckerPragerXXVE_UpdateDt( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context );
        void _DruckerPragerXXVE_Destroy( void* rheology, void* data );
        void _DruckerPragerXXVE_Build( void* rheology, void* data );
        void _DruckerPragerXXVE_Initialise( void* rheology, void* data );
        void _DruckerPragerXXVE_ModifyConstitutiveMatrix( 
				void*						rheology, 
				ConstitutiveMatrix*  		constitutiveMatrix,
				MaterialPointsSwarm* 		materialPointsSwarm,
				Element_LocalIndex			lElement_I,
				MaterialPoint*				materialPoint,
				Coord						xi 
		);

#endif

