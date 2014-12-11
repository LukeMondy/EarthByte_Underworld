#include <stdio.h>
#include <mpi.h>

#include "Components.h"

Bool Spherical_Components_Init( int* argc, char** argv[] ) {
   Stg_ComponentRegister* compReg = Stg_ComponentRegister_Get_ComponentRegister();
   Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

   Stg_ComponentRegister_Add( compReg, SphericalGenerator_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_SphericalGenerator_DefaultNew );
   RegisterParent( SphericalGenerator_Type, CartesianGenerator_Type );

   Stg_ComponentRegister_Add( compReg, RSGenerator_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_RSGenerator_DefaultNew );
   RegisterParent( RSGenerator_Type, CartesianGenerator_Type );

   Stg_ComponentRegister_Add( compReg, Q2SphericalGenerator_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)Q2SphericalGenerator_New );
   RegisterParent( Q2SphericalGenerator_Type, CartesianGenerator_Type );

   Stg_ComponentRegister_Add( compReg, ProjectionGenerator_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_ProjectionGenerator_DefaultNew );
   RegisterParent( ProjectionGenerator_Type, MeshGenerator_Type );

   Stg_ComponentRegister_Add( compReg, SphericalStiffnessMatrix_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)SphericalStiffnessMatrix_DefaultNew );
   RegisterParent( SphericalStiffnessMatrix_Type, StiffnessMatrix_Type );

   Stg_ComponentRegister_Add( compReg, SphericalForceVector_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_SphericalForceVector_DefaultNew );
   RegisterParent( SphericalForceVector_Type, ForceVector_Type );

   Stg_ComponentRegister_Add( compReg, Mesh_SphericalAlgorithms_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)Mesh_SphericalAlgorithms_New );
   RegisterParent( Mesh_SphericalAlgorithms_Type, Mesh_Algorithms_Type );

   Stg_ComponentRegister_Add( compReg, Mesh_ProjectionAlgorithms_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)Mesh_ProjectionAlgorithms_New );
   RegisterParent( Mesh_ProjectionAlgorithms_Type, Mesh_Algorithms_Type );

   Stg_ComponentRegister_Add( compReg, SphericalPeriodicAdvector_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_SphericalPeriodicAdvector_DefaultNew );
   RegisterParent( SphericalPeriodicAdvector_Type, PeriodicBoundariesManager_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_PointGravity_Type, (Name)"0", _Ppc_PointGravity_DefaultNew  );
   RegisterParent( Ppc_PointGravity_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Q2ProjectionGenerator_Type, (Name)"0", _Q2ProjectionGenerator_DefaultNew  );
   RegisterParent( Q2ProjectionGenerator_Type,CartesianGenerator_Type );

   Stg_ComponentRegister_Add( compReg, C0SphericalGenerator_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)C0SphericalGenerator_New );
   RegisterParent( C0SphericalGenerator_Type, C0Generator_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_ElementLine_Type, (Name)"0", _Ppc_ElementLine_DefaultNew  );
   RegisterParent( Ppc_ElementLine_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_Quality_Type, (Name)"0", _Ppc_Quality_DefaultNew  );
   RegisterParent( Ppc_Quality_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, LatLongRegion_Type, (Name)"0", _LatLongRegion_DefaultNew  );
   RegisterParent( LatLongRegion_Type, Stg_Shape_Type );

   Stg_ComponentRegister_Add( compReg, SphericalFeMesh_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_SphericalFeMesh_DefaultNew );
   RegisterParent( SphericalFeMesh_Type, FeMesh_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_Exponential_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_Exponential_DefaultNew );
   RegisterParent( Ppc_Exponential_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_2ndInvariant_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_2ndInvariant_DefaultNew );
   RegisterParent( Ppc_2ndInvariant_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_AdiabaticHeating_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_AdiabaticHeating_DefaultNew );
   RegisterParent( Ppc_AdiabaticHeating_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, MatAssembly_NA__Fi__NB_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_MatAssembly_NA__Fi__NB_DefaultNew );
   RegisterParent( MatAssembly_NA__Fi__NB_Type, StiffnessMatrixTerm_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_SphericalDepth_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_SphericalDepth_DefaultNew );
   RegisterParent( Ppc_SphericalDepth_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_Compression_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_Compression_DefaultNew );
   RegisterParent( Ppc_Compression_Type, Ppc_Type );

   Stg_ComponentRegister_Add( compReg, Ppc_VecDotVec_Type, (Name)"0", (Stg_Component_DefaultConstructorFunction*)_Ppc_VecDotVec_DefaultNew );
   RegisterParent( Ppc_VecDotVec_Type, Ppc_Type );

	VariableCondition_Register_Add( variableCondition_Register, SphereBC_Type, SphereBC_Factory );

   ConditionFunction_Register_Add( 
         condFunc_Register, 
         ConditionFunction_New( Cylindrical_WithPerturbation, (Name)"Cylindrical_WithPerturbation", NULL) 
         );

   return True;
}


