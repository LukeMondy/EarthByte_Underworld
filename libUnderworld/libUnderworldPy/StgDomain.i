/* -*- C -*-  (not really, but good for syntax highlighting) */

%module StgDomain

%import "StGermain.i"

%{
/* Includes the header in the wrapper code */
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
%}

#ifndef READ_HDF5
#define READ_HDF5
#endif

#ifndef WRITE_HDF5
#define WRITE_HDF5
#endif

%include "StgDomain_Typemaps.i"

%include "Geometry/ComplexMath.h"
%include "Geometry/ComplexVectorMath.h"
%include "Geometry/Dimension.h"
%include "Geometry/Edge.h"
%include "Geometry/Hex.h"
%include "Geometry/Line.h"
%include "Geometry/Plane.h"
%include "Geometry/QuadEdge.h"
%include "Geometry/RMatrix.h"
%include "Geometry/Simplex.h"
%include "Geometry/TensorMath.h"
%include "Geometry/TensorMultMath.h"
%include "Geometry/TrigMath.h"
%include "Geometry/types.h"
%include "Geometry/units.h"
%include "Geometry/VectorMath.h"
%include "Mesh/MeshClass.h"
%include "Mesh/MeshGenerator.h"
%include "Mesh/MeshAdaptor.h"
%include "Mesh/CartesianGenerator.h"
%include "Mesh/CompressionAdaptor.h"
%include "Mesh/GALESurfaceAdaptor.h"
%include "Mesh/MeshVariable.h"
%include "Mesh/Grid.h"
%include "Mesh/linearSpaceAdaptor.h"
%include "Mesh/Mesh_ElementType.h"
%include "Mesh/Mesh_Algorithms.h"
%include "Mesh/Mesh_CentroidAlgorithms.h"
%include "Mesh/Mesh_CentroidType.h"
%include "Mesh/Mesh_HexAlgorithms.h"
%include "Mesh/Mesh_HexType.h"
%include "Mesh/Mesh_RegularAlgorithms.h"
%include "Mesh/MeshTopology.h"
%include "Mesh/NodeBunching.h"
%include "Mesh/Remesher.h"
%include "Mesh/shortcuts.h"
%include "Mesh/SurfaceAdaptor.h"
%include "Mesh/types.h"
%include "Shape/ShapeClass.h"
%include "Shape/BelowPlane.h"
%include "Shape/BelowCosinePlane.h"
%include "Shape/Box.h"
%include "Shape/ConvexHull.h"
%include "Shape/Cylinder.h"
%include "Shape/Everywhere.h"
%include "Shape/Intersection.h"
%include "Shape/PolygonShape.h"
%include "Shape/PolygonShapeXZ.h"
%include "Shape/PythonShape.h"
%include "Shape/Sphere.h"
%include "Shape/Superellipsoid.h"
%include "Shape/types.h"
%include "Shape/Union.h"
%include "Swarm/ParticleLayout.h"
%include "Swarm/GlobalParticleLayout.h"
%include "Swarm/BackgroundParticleLayout.h"
%include "Swarm/SwarmVariable.h"
%include "Swarm/CellLayout.h"
%include "Swarm/ElementCellLayout.h"
%include "Swarm/CLLCellLayout.h"
%include "Swarm/LineParticleLayout.h"
%include "Swarm/ManualParticleLayout.h"
%include "Swarm/PerCellParticleLayout.h"
%include "Swarm/MeshParticleLayout.h"
%include "Swarm/OperatorSwarmVariable.h"
%include "Swarm/ParticleCommHandler.h"
%include "Swarm/ParticleMovementHandler.h"
%include "Swarm/ParticleShadowSync.h"
%include "Swarm/RandomParticleLayout.h"
%include "Swarm/ShadowInfo.h"
%include "Swarm/SingleCellLayout.h"
%include "Swarm/SpaceFillerParticleLayout.h"
%include "Swarm/StandardParticle.h"
%include "Swarm/Swarm_Register.h"
%include "Swarm/SwarmClass.h"
%include "Swarm/SwarmDump.h"
%include "Swarm/SwarmOutput.h"
%include "Swarm/SwarmVariable_Register.h"
%include "Swarm/TriGaussParticleLayout.h"
%include "Swarm/TriSingleCellLayout.h"
%include "Swarm/types.h"
%include "Swarm/UnionParticleLayout.h"
%include "Swarm/WithinShapeParticleLayout.h"
%include "Swarm/GaussParticleLayout.h"
%include "Swarm/GaussBorderParticleLayout.h"
%include "Swarm/IntegrationPoint.h"
%include "Swarm/PlaneParticleLayout.h"
%include "Utils/WallVC.h"
%include "Utils/AllElementsVC.h"
%include "Utils/AllNodesVC.h"
%include "Utils/ContactVC.h"
%include "Utils/FieldVariable.h"
%include "Utils/CoordinateFieldVariable.h"
%include "Utils/CornerVC.h"
%include "Utils/DofLayout.h"
%include "Utils/DomainContext.h"
%include "Utils/FieldVariable_Register.h"
%include "Utils/InnerWallVC.h"
%include "Utils/LinearRegression.h"
%include "Utils/MeshBoundaryShape.h"
%include "Utils/MeshShapeVC.h"
%include "Utils/NewRemesher.h"
%include "Utils/Operator.h"
%include "Utils/OperatorFieldVariable.h"
%include "Utils/RegRemesh.h"
%include "Utils/RegularMeshUtils.h"
%include "Utils/RegularRemesher.h"
%include "Utils/RegularRemesherCmpt.h"
%include "Utils/ShapeAdvector.h"
%include "Utils/SobolGenerator.h"
%include "Utils/TimeIntegrand.h"
%include "Utils/TimeIntegrator.h"
%include "Utils/types.h"
%include "Utils/Utils.h"

//
//%extend FieldVariable
//{
//
//  double* interResult=NULL;
//
//  int interpolateValueAt( Coord coord,  double* value )
//  {
//    return $self->_interpolateValueAt($self, coord, value);
//  }
//
//}
//