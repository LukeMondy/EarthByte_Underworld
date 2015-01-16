#ifndef __Spherical_Components_h__
#define __Spherical_Components_h__
	
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "rayCast.h"
#include "SphericalUtils.h"
#include "SphereBC.h"
#include "SphericalFeMesh.h"
#include "Mesh_SphericalAlgorithms.h"
#include "Mesh_ProjectionAlgorithms.h"
#include "SphericalGenerator.h"
#include "RSGenerator.h"
#include "ProjectionGenerator.h"
#include "C0SphericalGenerator.h"
#include "SphericalStiffnessMatrix.h"
#include "SphericalForceVector.h"
#include "SphericalPeriodicAdvector.h"
#include "Q2SphericalGenerator.h"
#include "Q2ProjectionGenerator.h"
#include "Ppc_PointGravity.h"
#include "Ppc_ElementLine.h"
#include "Ppc_Quality.h"
#include "LatLongRegion.h"
#include "Ppc_Exponential.h"
#include "Ppc_AdiabaticHeating.h"
#include "MatAssembly_NA__Fi__NB.h"
#include "Ppc_SphericalDepth.h"
#include "Ppc_TensorInvariant.h"
#include "Ppc_VecDotVec.h"
#include "SLIntegrator_Spherical.h"
#include "SLIntegrator_Polar.h"
#include "SLIntegrator_FullSphere.h"

#include "Init.h"
#include "Finalise.h"

#endif 
