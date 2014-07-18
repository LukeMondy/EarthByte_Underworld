/* -*- C -*-  (not really, but good for syntax highlighting) */

%module ImportersToolbox

%import "StGermain.i"
%import "StgDomain.i"
%import "StgFEM.i"
%import "PICellerator.i"
%import "Underworld.i"
%import "gLucifer.i"

%{
/* Includes the header in the wrapper code */
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <ImportersToolbox/ImportersToolbox.h>
#include <gLucifer/gLucifer.h>
%}

#%include "Base/VoxelParticleLayout.h"
#%include "Base/VoxelHeightFieldVariable.h"
#%include "Base/VoxelFieldVariable_GMT.h"
#%include "Base/VoxelDataHandler_VTKStructuredPoints.h"
#%include "Base/VoxelDataHandler_GocadProperties.h"
#%include "Base/VoxelDataHandler_GocadMaterials.h"
#%include "Base/VoxelDataHandler_GocadAbstract.h"
#%include "Base/VoxelDataHandler_Geomodeller.h"
#%include "Base/VoxelDataHandler_GMT.h"
#%include "Base/VoxelDataHandler_Abstract.h"
#%include "Base/SpatialDataFieldVariable.h"
%include "Base/VoxelFieldVariable.h"
#%include "Base/VoxelDataHandler_ASCII.h"
#%include "Base/VoxelDataHandler_ndarray.h"
#%include "Base/types.h"

