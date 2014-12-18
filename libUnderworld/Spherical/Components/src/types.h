#ifndef __Spherical_Components_types_h__
#define __Spherical_Components_types_h__

typedef struct Mesh_SphericalAlgorithms Mesh_SphericalAlgorithms;
typedef struct Mesh_ProjectionAlgorithms Mesh_ProjectionAlgorithms;
typedef struct SphericalGenerator SphericalGenerator;
typedef struct RSGenerator RSGenerator;
typedef struct Q2SphericalGenerator Q2SphericalGenerator;
typedef struct Q2ProjectionGenerator Q2ProjectionGenerator;
typedef struct ProjectionGenerator ProjectionGenerator;
typedef struct C0SphericalGenerator C0SphericalGenerator;
typedef struct SphericalStiffnessMatrix SphericalStiffnessMatrix;
typedef struct SphericalForceVector SphericalForceVector;
typedef struct SphericalPeriodicAdvector SphericalPeriodicAdvector;
typedef struct Ppc_PointGravity Ppc_PointGravity;
typedef struct Ppc_ElementLine Ppc_ElementLine;
typedef struct Ppc_Quality Ppc_Quality;
typedef struct LatLongRegion LatLongRegion;
typedef struct SphericalFeMesh SphericalFeMesh;
typedef struct SphereBC SphereBC;
typedef struct Ppc_Exponential Ppc_Exponential;
typedef struct Ppc_AdiabaticHeating Ppc_AdiabaticHeating;
typedef struct MatAssembly_NA__Fi__NB MatAssembly_NA__Fi__NB;
typedef struct Ppc_SphericalDepth Ppc_SphericalDepth;
typedef struct Ppc_TensorInvariant Ppc_TensorInvariant;
typedef struct Ppc_VecDotVec Ppc_VecDotVec;
typedef struct SLIntegrator_Spherical SLIntegrator_Spherical;

typedef enum {
   SphereBC_Wall_Inner,
   SphereBC_Wall_Outer,
   SphereBC_Wall_Size
} SphereBC_Wall;


#endif /* __Spherical_Components_types_h__ */
