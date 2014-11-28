#ifndef __Spherical_Components_SphericalStiffnessMatrix_h__
#define __Spherical_Components_SphericalStiffnessMatrix_h__
	

	/* Textual name of this class */
	extern const Type SphericalStiffnessMatrix_Type;
	
	/* StiffnessMatrix information */
	#define __SphericalStiffnessMatrix  \
		/* General info */ \
		__StiffnessMatrix \
									\
		double* rotMat; /* storage for rotation matrix */ \
		double* tmpMat; /* storage for rotation matrix calculation, DGEMM needs a third matrix*/ \

	struct SphericalStiffnessMatrix { __SphericalStiffnessMatrix };
	
	/* Creation implementation / Virtual constructor */
	void* SphericalStiffnessMatrix_DefaultNew( Name name );
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SPHERICALSTIFFNESSMATRIX_DEFARGS \
                STIFFNESSMATRIX_DEFARGS

	#define SPHERICALSTIFFNESSMATRIX_PASSARGS \
                STIFFNESSMATRIX_PASSARGS

	SphericalStiffnessMatrix* _SphericalStiffnessMatrix_New(  STIFFNESSMATRIX_DEFARGS  );		
		
	void _SphericalStiffnessMatrix_Init( 
			SphericalStiffnessMatrix*                                 self);
	
	void _SphericalStiffnessMatrix_AssignFromXML( void* stiffnessMatrix, Stg_ComponentFactory* cf, void* data );

	/* Initialisation implementation */
	void _SphericalStiffnessMatrix_Build( void* stiffnessMatrix, void* data );
	
	/* Destruction implementation */
	void _SphericalStiffnessMatrix_Destroy( void* stiffnessMatrix, void* data );

	/* important function. After matrix assemble it will apply rotation operations to properly specify non axis alined boundary condition dofs */
	void SphericalStiffnessMatrix_AssembleElement(
		void* stiffnessMatrix,
		Element_LocalIndex element_lI,
		SystemLinearEquations* sle,
		FiniteElementContext* context,
		double** elStiffMatVecToAdd);


	/* TODO: don't think this function is really needed, __StiffnessMatrix_NewAssemble, should be fine. But for development I leave this in */
	void _SphericalStiffnessMatrix_NewAssemble( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context );

  /** calculate the rotation matrix for this element, assume it is already on the boundary of the domain **/
	Bool SphericalStiffnessMatrix_EvaluateRotationMatrix2D( FeVariable* var, unsigned element_lI, double *rMat );
	Bool SphericalStiffnessMatrix_EvaluateRotationMatrix3D( FeVariable* var, unsigned element_lI, double *rMat );
	
#endif 
