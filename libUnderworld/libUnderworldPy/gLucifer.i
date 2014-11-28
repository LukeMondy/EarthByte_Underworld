/* -*- C -*-  (not really, but good for syntax highlighting) */

%module gLucifer

%import "StGermain.i"
%import "StgDomain.i"
%import "StgFEM.i"
%import "PICellerator.i"
%import "Underworld.i"

%{
/* Includes the header in the wrapper code */
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <gLucifer/gLucifer.h>
%}

%include "Base/Base.h"
%include "Base/Window.h"
%include "Base/Camera.h"
%include "Base/ColourMap.h"
%include "Base/Database.h"
%include "Base/DrawingObject.h"
%include "Base/DrawingObject_Register.h"
%include "Base/Finalise.h"
%include "Base/Init.h"
%include "Base/sqlite3.h"
%include "Base/types.h"
%include "Base/Viewport.h"
%include "Base/X11Colours.h"
%include "DrawingObjects/SwarmViewer.h"
%include "DrawingObjects/Axis.h"
%include "DrawingObjects/Capture.h"
%include "DrawingObjects/ColourBar.h"
%include "DrawingObjects/Contour.h"
%include "DrawingObjects/ContourCrossSection.h"
%include "DrawingObjects/CrossSection.h"
%include "DrawingObjects/DrawingObjects.h"
%include "DrawingObjects/Eigenvectors.h"
%include "DrawingObjects/EigenvectorsCrossSection.h"
%include "DrawingObjects/FeVariableSurface.h"
%include "DrawingObjects/FieldSampler.h"
%include "DrawingObjects/FieldVariableBorder.h"
%include "DrawingObjects/Finalise.h"
%include "DrawingObjects/HistoricalSwarmTrajectory.h"
%include "DrawingObjects/Init.h"
%include "DrawingObjects/Isosurface.h"
%include "DrawingObjects/IsosurfaceCrossSection.h"
%include "DrawingObjects/MeshCrossSection.h"
%include "DrawingObjects/MeshViewer.h"
%include "DrawingObjects/Plot.h"
%include "DrawingObjects/ScalarField.h"
%include "DrawingObjects/ScalarFieldCrossSection.h"
%include "DrawingObjects/ScalarFieldOnMesh.h"
%include "DrawingObjects/ScalarFieldOnMeshCrossSection.h"
%include "DrawingObjects/SwarmRGBColourViewer.h"
%include "DrawingObjects/SwarmShapes.h"
%include "DrawingObjects/SwarmVectors.h"        
%include "DrawingObjects/TextureMap.h"
%include "DrawingObjects/TimeStep.h"
%include "DrawingObjects/Title.h"
%include "DrawingObjects/types.h"
%include "DrawingObjects/VectorArrowCrossSection.h"
%include "DrawingObjects/VectorArrowMeshCrossSection.h"
%include "DrawingObjects/VectorArrows.h"
%include "DrawingObjects/VectorArrowsOnMesh.h"
%include "libgLucifer/Finalise.h"
%include "libgLucifer/gLucifer.h"
%include "libgLucifer/Init.h"
%include "libgLucifer/Toolbox/Toolbox.h"



