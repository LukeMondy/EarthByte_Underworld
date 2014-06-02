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

#ifndef __Experimental_Rheology_ViscousPenaltyConstMatrixCartesian_h__
#define __Experimental_Rheology_ViscousPenaltyConstMatrixCartesian_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type ViscousPenaltyConstMatrixCartesian_Type;

	/** ViscousPenaltyConstMatrixCartesian class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __ViscousPenaltyConstMatrixCartesian \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ConstitutiveMatrix \
		\
		/* Virtual functions go here */ \
		\
		/* ViscousPenaltyConstMatrixCartesian info */ \
		double** Dtilda_B; \
		double*  Ni; \
		/* Parameter */ \
		double   incompressibility_Penalty; \
	    Bool     viscosityWeighting;

	struct ViscousPenaltyConstMatrixCartesian { __ViscousPenaltyConstMatrixCartesian };

	ViscousPenaltyConstMatrixCartesian* ViscousPenaltyConstMatrixCartesian_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffnessMatrix,
		Swarm*                                              swarm,
		Dimension_Index                                     dim,
		FiniteElementContext*                               context,	
		Materials_Register*                                 materials_Register,
		double                                              incompressibility_Penalty,
		Bool                                                viscosityWeighting
		);

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define VISCOUSPENALTYCONSTMATRIXCARTESIAN_DEFARGS \
                CONSTITUTIVEMATRIX_DEFARGS

	#define VISCOUSPENALTYCONSTMATRIXCARTESIAN_PASSARGS \
                CONSTITUTIVEMATRIX_PASSARGS

	ViscousPenaltyConstMatrixCartesian* _ViscousPenaltyConstMatrixCartesian_New(  VISCOUSPENALTYCONSTMATRIXCARTESIAN_DEFARGS  );
	
   void _ViscousPenaltyConstMatrixCartesian_Init( 
		ViscousPenaltyConstMatrixCartesian*         self,
		double                                      incompressibility_Penalty,
		Bool										viscosityWeighting );
	
	void  _ViscousPenaltyConstMatrixCartesian_Delete( void* constitutiveMatrix );
	void  _ViscousPenaltyConstMatrixCartesian_Print( void* constitutiveMatrix, Stream* stream );

	void* _ViscousPenaltyConstMatrixCartesian_DefaultNew( Name name ) ;
	void  _ViscousPenaltyConstMatrixCartesian_AssembleFromXML( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) ;
	void  _ViscousPenaltyConstMatrixCartesian_Build( void* constitutiveMatrix, void* data ) ;
	void  _ViscousPenaltyConstMatrixCartesian_Initialise( void* constitutiveMatrix, void* data ) ;
	void  _ViscousPenaltyConstMatrixCartesian_Execute( void* constitutiveMatrix, void* data ) ;
	void  _ViscousPenaltyConstMatrixCartesian_Destroy( void* constitutiveMatrix, void* data ) ;

	void  _ViscousPenaltyConstMatrixCartesian_AssembleElement( 
		void*                                                constitutiveMatrix,
		StiffnessMatrix*                                     stiffnessMatrix, 
		Element_LocalIndex                                   lElement_I, 
		SystemLinearEquations*                               sle,
		FiniteElementContext*                                context,
		double**                                             elStiffMat ) ;

	void _ViscousPenaltyConstMatrixCartesian2D_SetValueInAllEntries( void* constitutiveMatrix, double value ) ;
	void _ViscousPenaltyConstMatrixCartesian3D_SetValueInAllEntries( void* constitutiveMatrix, double value ) ;

	double _ViscousPenaltyConstMatrixCartesian2D_GetIsotropicViscosity( void* constitutiveMatrix ) ;
	double _ViscousPenaltyConstMatrixCartesian3D_GetIsotropicViscosity( void* constitutiveMatrix ) ;

	void _ViscousPenaltyConstMatrixCartesian2D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) ;
	void _ViscousPenaltyConstMatrixCartesian3D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) ;

	void _ViscousPenaltyConstMatrixCartesian2D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director );
	void _ViscousPenaltyConstMatrixCartesian3D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director );

	void _ViscousPenaltyConstMatrixCartesian2D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B );
	void _ViscousPenaltyConstMatrixCartesian3D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B );

	void _ViscousPenaltyConstMatrixCartesian2D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) ;
	void _ViscousPenaltyConstMatrixCartesian3D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) ;

	/* a function which defines the storage of each particle's constitutive information on the particle, 
	 * should be called before the "Build" phase */
	void ViscousPenaltyConstMatrixCartesian_SetupParticleStorage( ViscousPenaltyConstMatrixCartesian* self );

#endif

