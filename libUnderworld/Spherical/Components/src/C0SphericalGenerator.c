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


/* Textual name of this class */
const Type C0SphericalGenerator_Type = "C0SphericalGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

C0SphericalGenerator* C0SphericalGenerator_New( Name name, AbstractContext* context ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(C0SphericalGenerator);
	Type                                                      type = C0SphericalGenerator_Type;
	Stg_Class_DeleteFunction*                              _delete = _C0Generator_Delete;
	Stg_Class_PrintFunction*                                _print = _C0Generator_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (void* (*)(Name))_C0Generator_New;
	Stg_Component_ConstructFunction*                    _construct = _C0Generator_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _C0Generator_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _C0Generator_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _C0Generator_Execute;
	Stg_Component_DestroyFunction*                        _destroy = NULL;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	MeshGenerator_SetDimSizeFunc*                   setDimSizeFunc = _MeshGenerator_SetDimSize;
	MeshGenerator_GenerateFunc*                       generateFunc = (MeshGenerator_GenerateFunc*)C0SphericalGenerator_Generate;

	C0SphericalGenerator* self = (C0SphericalGenerator*)_C0Generator_New(  C0GENERATOR_PASSARGS  );
   
   _MeshGenerator_Init( (MeshGenerator*)self, context );
	_C0Generator_Init( (C0Generator*)self );

   return self;
}

void C0SphericalGenerator_Generate( void* generator, void* _mesh ) {
	C0SphericalGenerator*	self = (C0SphericalGenerator*)generator;
	FeMesh*		mesh = (FeMesh*)_mesh;

	assert( self && Stg_CheckType( self, C0SphericalGenerator ) );
	assert( mesh && Stg_CheckType( mesh, FeMesh ) );

	C0Generator_BuildTopology( (C0Generator*)self, mesh );
	C0Generator_BuildGeometry( (C0Generator*) self, mesh );
	C0Generator_BuildElementTypes( (C0Generator*) self, mesh );

	mesh->elementMesh = True; // ie this mesh defines the elements, not the nodes
}

