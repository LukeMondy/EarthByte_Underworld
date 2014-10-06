/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"
#include "Q2ProjectionGenerator.h"


/* Textual name of this class */
const Type Q2ProjectionGenerator_Type = "Q2ProjectionGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

void* _Q2ProjectionGenerator_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                    _sizeOfSelf = sizeof(Q2ProjectionGenerator);
	Type                                                            type = Q2ProjectionGenerator_Type;
	Stg_Class_DeleteFunction*                                    _delete = _C2Generator_Delete;
	Stg_Class_PrintFunction*                                      _print = _C2Generator_Print;
	Stg_Class_CopyFunction*                                        _copy = NULL;
	Stg_Component_DefaultConstructorFunction*        _defaultConstructor = (void* (*)(Name))_Q2ProjectionGenerator_New;
	Stg_Component_ConstructFunction*                          _construct = _ProjectionGenerator_AssignFromXML;
	Stg_Component_BuildFunction*                                  _build = _Q2ProjectionGenerator_Build;
	Stg_Component_InitialiseFunction*                        _initialise = _ProjectionGenerator_Initialise;
	Stg_Component_ExecuteFunction*                              _execute = _C2Generator_Execute;
	Stg_Component_DestroyFunction*                              _destroy = NULL;
	AllocationType                                    nameAllocationType = NON_GLOBAL;
	MeshGenerator_SetDimSizeFunc*                         setDimSizeFunc = CartesianGenerator_SetDimSize;
	MeshGenerator_GenerateFunc*                             generateFunc = CartesianGenerator_Generate;
	CartesianGenerator_SetTopologyParamsFunc*      setTopologyParamsFunc = C2Generator_SetTopologyParams;
	CartesianGenerator_GenElementsFunc*                  genElementsFunc = _CartesianGenerator_GenElements;
	CartesianGenerator_GenFacesFunc*                        genFacesFunc = _CartesianGenerator_GenFaces;
	CartesianGenerator_GenEdgesFunc*                        genEdgesFunc = _CartesianGenerator_GenEdges;
	CartesianGenerator_GenVerticesFunc*                  genVerticesFunc = _CartesianGenerator_GenVertices;
	CartesianGenerator_GenElementVertexIncFunc*  genElementVertexIncFunc = C2Generator_GenElementVertexInc;
	CartesianGenerator_GenVolumeEdgeIncFunc*        genVolumeEdgeIncFunc = _CartesianGenerator_GenVolumeEdgeInc;
	CartesianGenerator_GenVolumeFaceIncFunc*        genVolumeFaceIncFunc = _CartesianGenerator_GenVolumeFaceInc;
	CartesianGenerator_GenFaceVertexIncFunc*        genFaceVertexIncFunc = C2Generator_GenFaceVertexInc;
	CartesianGenerator_GenFaceEdgeIncFunc*            genFaceEdgeIncFunc = _CartesianGenerator_GenFaceEdgeInc;
	CartesianGenerator_GenEdgeVertexIncFunc*        genEdgeVertexIncFunc = C2Generator_GenEdgeVertexInc;
	CartesianGenerator_GenElementTypesFunc*          genElementTypesFunc = C2Generator_GenElementTypes;
	CartesianGenerator_CalcGeomFunc*                calcGeomFunc = NULL;

	Q2ProjectionGenerator* self = _Q2ProjectionGenerator_New(  C2SPHERICALGENERATOR_PASSARGS  );

   return self;
}

Q2ProjectionGenerator* _Q2ProjectionGenerator_New(  C2SPHERICALGENERATOR_DEFARGS  ) {
	Q2ProjectionGenerator*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Q2ProjectionGenerator) );
	self = (Q2ProjectionGenerator*)_CartesianGenerator_New(  CARTESIANGENERATOR_PASSARGS  );

	return self;
}

void _Q2ProjectionGenerator_Delete( void* generator ) {
	_C2Generator_Delete( generator );
}

void _Q2ProjectionGenerator_Build( void* meshGenerator, void* data ) {
	Q2ProjectionGenerator*		self = (Q2ProjectionGenerator*)meshGenerator;

	_ProjectionGenerator_Build( meshGenerator, data );

	// we are creating a quadratic node layout so we increase the number of 'squares' in the algorithm
	// TODO:check in 3D
	self->nSquares = 2 * self->nSquares - 1;
}
