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

#include "Components.h"


/* Textual name of this class */
const Type Mesh_SphericalAlgorithms_Type = "Mesh_SphericalAlgorithms";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_SphericalAlgorithms* Mesh_SphericalAlgorithms_New( Name name, AbstractContext* context ) {
	/* Variables set in this function */
	SizeT                                                   _sizeOfSelf = sizeof(Mesh_SphericalAlgorithms);
	Type                                                           type = Mesh_SphericalAlgorithms_Type;
	Stg_Class_DeleteFunction*                                   _delete = _Mesh_SphericalAlgorithms_Delete;
	Stg_Class_PrintFunction*                                     _print = _Mesh_SphericalAlgorithms_Print;
	Stg_Class_CopyFunction*                                       _copy = NULL;
	Stg_Component_DefaultConstructorFunction*       _defaultConstructor = (void* (*)(Name))_Mesh_SphericalAlgorithms_New;
	Stg_Component_ConstructFunction*                         _construct = _Mesh_SphericalAlgorithms_AssignFromXML;
	Stg_Component_BuildFunction*                                 _build = _Mesh_SphericalAlgorithms_Build;
	//Stg_Component_InitialiseFunction*                       _initialise = _Mesh_SphericalAlgorithms_Initialise;
	Stg_Component_InitialiseFunction*                       _initialise = _Mesh_Algorithms_Initialise;
	Stg_Component_ExecuteFunction*                             _execute = _Mesh_SphericalAlgorithms_Execute;
	Stg_Component_DestroyFunction*                             _destroy = _Mesh_SphericalAlgorithms_Destroy;
	AllocationType                                   nameAllocationType = NON_GLOBAL;
	Mesh_Algorithms_SetMeshFunc*                            setMeshFunc = Mesh_SphericalAlgorithms_SetMesh;
	Mesh_Algorithms_UpdateFunc*                              updateFunc = Mesh_SphericalAlgorithms_Update;
	Mesh_Algorithms_NearestVertexFunc*                nearestVertexFunc = _Mesh_Algorithms_NearestVertex;
	Mesh_Algorithms_SearchFunc*                              searchFunc = _Mesh_Algorithms_Search;
	Mesh_Algorithms_SearchElementsFunc*              searchElementsFunc = _Mesh_Algorithms_SearchElements;
//	Mesh_Algorithms_SearchElementsFunc*              searchElementsFunc = Mesh_SphericalAlgorithms_SearchElements; 
	Mesh_Algorithms_GetMinimumSeparationFunc*  getMinimumSeparationFunc = _Mesh_Algorithms_GetMinimumSeparation;
	Mesh_Algorithms_GetLocalCoordRangeFunc*      getLocalCoordRangeFunc = _Mesh_SphericalAlgorithms_GetLocalCoordRange;
	Mesh_Algorithms_GetDomainCoordRangeFunc*    getDomainCoordRangeFunc = _Mesh_Algorithms_GetDomainCoordRange;
	Mesh_Algorithms_GetGlobalCoordRangeFunc*    getGlobalCoordRangeFunc = _Mesh_SphericalAlgorithms_GetGlobalCoordRange;

	Mesh_SphericalAlgorithms* self = _Mesh_SphericalAlgorithms_New(  MESH_REGULARALGORITHMS_PASSARGS  );

	/* Mesh_SphericalAlgorithms info */
	_Mesh_Algorithms_Init( (Mesh_Algorithms*)self, context );
	_Mesh_SphericalAlgorithms_Init( self );

   return self;
}

Mesh_SphericalAlgorithms* _Mesh_SphericalAlgorithms_New(  MESH_REGULARALGORITHMS_DEFARGS  ) {
	Mesh_SphericalAlgorithms* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Mesh_SphericalAlgorithms) );
	self = (Mesh_SphericalAlgorithms*)_Mesh_Algorithms_New(  MESH_ALGORITHMS_PASSARGS  );

	return self;
}

void _Mesh_SphericalAlgorithms_Init( void* algorithms ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_SphericalAlgorithms ) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_SphericalAlgorithms_Delete( void* algorithms ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _Mesh_SphericalAlgorithms_Print( void* algorithms, Stream* stream ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;
	
	/* Print parent */
	Journal_Printf( stream, "Mesh_SphericalAlgorithms (ptr): (%p)\n", self );
	_Mesh_Algorithms_Print( self, stream );
}

void _Mesh_SphericalAlgorithms_AssignFromXML( void* algorithms, Stg_ComponentFactory* cf, void* data ) {
	_Mesh_Algorithms_AssignFromXML( algorithms, cf, data );
   _Mesh_SphericalAlgorithms_Init( algorithms );
}

void _Mesh_SphericalAlgorithms_Build( void* algorithms, void* data ) {
    _Mesh_Algorithms_Build( algorithms, data );
}

void _Mesh_SphericalAlgorithms_Initialise( void* algorithms, void* data ) {
	Mesh_SphericalAlgorithms* self = (Mesh_SphericalAlgorithms*)algorithms;

	SphericalGenerator* sgen = NULL; 
	if( Stg_Class_IsInstance( self->mesh->generator, SphericalGenerator_Type ) ) {
		sgen = (SphericalGenerator*)self->mesh->generator;
	} else if( Stg_Class_IsInstance( self->mesh->generator, Q2SphericalGenerator_Type ) ) {
		sgen = (Q2SphericalGenerator*)self->mesh->generator;
	} else {
		assert(0); //error
	}

	_Mesh_Algorithms_Initialise( algorithms, data );

	// copy over spherical mesh infomation
	memcpy( self->sres, sgen->sph_res, 3*sizeof(double) );
	memcpy( self->smincrd, sgen->crdMin, 3*sizeof(double) );
	memcpy( self->smaxcrd, sgen->crdMax, 3*sizeof(double) );
}

void _Mesh_SphericalAlgorithms_Execute( void* algorithms, void* data ) {
    _Mesh_Algorithms_Execute( algorithms, data );
}

void _Mesh_SphericalAlgorithms_Destroy( void* algorithms, void* data ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;

	Mesh_SphericalAlgorithms_Destruct( self );

   _Mesh_Algorithms_Destroy( algorithms, data );
}

void Mesh_SphericalAlgorithms_SetMesh( void* algorithms, void* mesh ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_SphericalAlgorithms ) );

	Mesh_SphericalAlgorithms_Destruct( self );
	_Mesh_Algorithms_SetMesh( self, mesh );
}

void Mesh_SphericalAlgorithms_Update( void* algorithms ) {
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_SphericalAlgorithms ) );
	assert( self->mesh );

	Mesh_SphericalAlgorithms_Destruct( self );
	_Mesh_Algorithms_Update( self );

}

void _Mesh_SphericalAlgorithms_GetLocalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nVerts, nEls;
	unsigned*		verts;
	double*			vert, rtp[3];
	unsigned		nDims;
	unsigned		v_i, e_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	nEls = Mesh_GetLocalSize( mesh, MT_VERTEX );

  // get the first node as an estimate
  vert = Mesh_GetVertex( mesh, 0 );
  Spherical_XYZ2RTP2D( vert, rtp ); Spherical_Radians2Output( rtp );
	memcpy( min, rtp, nDims * sizeof(double) );
	memcpy( max, rtp, nDims * sizeof(double) );

	// loop over each node
	for( v_i = 1; v_i < nEls; v_i++ ) {
		vert = Mesh_GetVertex( mesh, v_i );
		Spherical_XYZ2RTP2D( vert, rtp ); Spherical_Radians2Output( rtp );
		// test vertex coords
		for( d_i = 0; d_i < nDims; d_i++ ) {
			if( rtp[d_i] < min[d_i] )
				min[d_i] = rtp[d_i];
			if( rtp[d_i] > max[d_i] )
				max[d_i] = rtp[d_i];
		}
	}
}

#if 0

void _Mesh_SphericalAlgorithms_GetDomainCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nVerts;
	double*			vert;
	unsigned		nDims;
	unsigned		v_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	nVerts = Mesh_GetDomainSize( mesh, MT_VERTEX );
	memcpy( min, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	memcpy( max, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	for( v_i = 1; v_i < nVerts; v_i++ ) {
		vert = Mesh_GetVertex( mesh, v_i );
		Mesh_SphericalAlgorithms_XYZ2RTP( vert );
		for( d_i = 0; d_i < nDims; d_i++ ) {
			if( vert[d_i] < min[d_i] )
				min[d_i] = vert[d_i];
			if( vert[d_i] > max[d_i] )
				max[d_i] = vert[d_i];
		}
	}
}

#endif
void _Mesh_SphericalAlgorithms_GetGlobalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nDims;
	double			*localMin, *localMax;
	MPI_Comm		comm;
	unsigned		d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	localMin = Memory_Alloc_Array_Unnamed( double, nDims );
	localMax = Memory_Alloc_Array_Unnamed( double, nDims );

	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
	Mesh_Algorithms_GetLocalCoordRange( self, localMin, localMax );
	for( d_i = 0; d_i < Mesh_GetDimSize( mesh ); d_i++ ) {
		MPI_Allreduce( localMin + d_i, min + d_i, 1, MPI_DOUBLE, MPI_MIN, comm );
		MPI_Allreduce( localMax + d_i, max + d_i, 1, MPI_DOUBLE, MPI_MAX, comm );
	}

	FreeArray( localMin );
	FreeArray( localMax );
}


Bool _Mesh_SphericalAlgorithms_Search( void* algorithms, void* _mesh, double* point, 
				     MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_SphericalAlgorithms*	self = (Mesh_SphericalAlgorithms*)algorithms;
	Mesh*			mesh = (Mesh*)_mesh;

	assert( self );
	assert( mesh );
	assert( dim );
	assert( ind );

	/* TODO */
	abort();

	return False;
}

Bool Mesh_SphericalAlgorithms_SearchElements( void* algorithms, double* xyz, unsigned* elInd ) {
   /*
      This algorithm is too precise for linear elements that don't actually map into a sperical segment.
      Particles near the edge of the inner and outer (radial) surfaces are mistakenly found to be outside or inside the element
   */
	Mesh_SphericalAlgorithms *self = (Mesh_SphericalAlgorithms*)algorithms;
	Mesh                     *mesh=NULL;
        SphericalGenerator       *sgen=NULL;
	unsigned                 nDims;
	unsigned		inds[3];
	Grid			*elGrid;
	double			out, frac, integer;
	unsigned		d_i;
	double r, theta, tmp, dim=2;
	double point[3];

	assert( self );
	assert( elInd );
	assert( Mesh_GetDimSize( self->mesh ) <= 3 );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	elGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, mesh->elGridId );
        sgen = (SphericalGenerator*)mesh->generator;

	/* translate point from cartesian coords, xyz, into polar coord, point 
		 r = sqrt( x^2 + y+2 + z^2 )
		 theta = atan( y/x )
	 */

	// convert xyz coords to rtp in degrees
        if( nDims == 2 )
           Spherical_XYZ2RTP2D( xyz, point );
        else
           Spherical_XYZ2RTP3D( xyz, point );

        // convert to radians
        Spherical_Radians2Output( point );

	for( d_i = 0; d_i < nDims; d_i++  ) {
		if( Num_Approx( point[d_i] - sgen->crdMax[d_i], 0.0 ) )
			inds[d_i] = elGrid->sizes[d_i] - 1;
		else if( point[d_i] < sgen->crdMin[d_i] || point[d_i] > sgen->crdMax[d_i] )
			return False;
		else {
			out = (point[d_i] - sgen->crdMin[d_i]) / sgen->sph_res[d_i];
			frac = modf( out, &integer );
			inds[d_i] = (unsigned)integer;
			if( inds[d_i] > 0 && Num_Approx( frac, 0.0 ) )
				inds[d_i]--;
		}
	}

	*elInd = Grid_Project( elGrid, inds );
	return Mesh_GlobalToDomain( mesh, nDims, *elInd, elInd );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Mesh_SphericalAlgorithms_Destruct( Mesh_SphericalAlgorithms* self ) {
	assert( self && Stg_CheckType( self, Mesh_SphericalAlgorithms ) );
}


