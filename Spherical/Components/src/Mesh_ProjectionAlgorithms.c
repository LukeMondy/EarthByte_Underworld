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
#include <float.h>
#include <mpi.h>


#include "Components.h"


/* Textual name of this class */
const Type Mesh_ProjectionAlgorithms_Type = "Mesh_ProjectionAlgorithms";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_ProjectionAlgorithms* Mesh_ProjectionAlgorithms_New( Name name, AbstractContext* context ) {
	/* Variables set in this function */
	SizeT                                                   _sizeOfSelf = sizeof(Mesh_ProjectionAlgorithms);
	Type                                                           type = Mesh_ProjectionAlgorithms_Type;
	Stg_Class_DeleteFunction*                                   _delete = _Mesh_ProjectionAlgorithms_Delete;
	Stg_Class_PrintFunction*                                     _print = _Mesh_ProjectionAlgorithms_Print;
	Stg_Class_CopyFunction*                                       _copy = NULL;
	Stg_Component_DefaultConstructorFunction*       _defaultConstructor = (void* (*)(Name))_Mesh_ProjectionAlgorithms_New;
	Stg_Component_ConstructFunction*                         _construct = _Mesh_ProjectionAlgorithms_AssignFromXML;
	Stg_Component_BuildFunction*                                 _build = _Mesh_ProjectionAlgorithms_Build;
	Stg_Component_InitialiseFunction*                       _initialise = _Mesh_ProjectionAlgorithms_Initialise;
	//Stg_Component_InitialiseFunction*                       _initialise = _Mesh_Algorithms_Initialise;
	Stg_Component_ExecuteFunction*                             _execute = _Mesh_ProjectionAlgorithms_Execute;
	Stg_Component_DestroyFunction*                             _destroy = _Mesh_ProjectionAlgorithms_Destroy;
	AllocationType                                   nameAllocationType = NON_GLOBAL;
	Mesh_Algorithms_SetMeshFunc*                            setMeshFunc = Mesh_ProjectionAlgorithms_SetMesh;
	Mesh_Algorithms_UpdateFunc*                              updateFunc = Mesh_ProjectionAlgorithms_Update;
	Mesh_Algorithms_NearestVertexFunc*                nearestVertexFunc = _Mesh_Algorithms_NearestVertex;
	Mesh_Algorithms_SearchFunc*                              searchFunc = _Mesh_Algorithms_Search;
	Mesh_Algorithms_SearchElementsFunc*              searchElementsFunc = _Mesh_Algorithms_SearchElements;
	//Mesh_Algorithms_SearchElementsFunc*              searchElementsFunc = Mesh_ProjectionAlgorithms_SearchElements;
	Mesh_Algorithms_GetMinimumSeparationFunc*  getMinimumSeparationFunc = _Mesh_Algorithms_GetMinimumSeparation;
	Mesh_Algorithms_GetLocalCoordRangeFunc*      getLocalCoordRangeFunc = _Mesh_Algorithms_GetLocalCoordRange;
	Mesh_Algorithms_GetDomainCoordRangeFunc*    getDomainCoordRangeFunc = _Mesh_Algorithms_GetDomainCoordRange;
	Mesh_Algorithms_GetGlobalCoordRangeFunc*    getGlobalCoordRangeFunc = _Mesh_Algorithms_GetGlobalCoordRange;

	Mesh_ProjectionAlgorithms* self = _Mesh_ProjectionAlgorithms_New(  MESH_REGULARALGORITHMS_PASSARGS  );

	/* Mesh_ProjectionAlgorithms info */
	_Mesh_Algorithms_Init( (Mesh_Algorithms*)self, context );
	_Mesh_ProjectionAlgorithms_Init( self );

   return self;
}

Mesh_ProjectionAlgorithms* _Mesh_ProjectionAlgorithms_New(  MESH_REGULARALGORITHMS_DEFARGS  ) {
	Mesh_ProjectionAlgorithms* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Mesh_ProjectionAlgorithms) );
	self = (Mesh_ProjectionAlgorithms*)_Mesh_Algorithms_New(  MESH_ALGORITHMS_PASSARGS  );

	return self;
}

void _Mesh_ProjectionAlgorithms_Init( void* algorithms ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_ProjectionAlgorithms ) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_ProjectionAlgorithms_Delete( void* algorithms ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _Mesh_ProjectionAlgorithms_Print( void* algorithms, Stream* stream ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;
	
	/* Print parent */
	Journal_Printf( stream, "Mesh_ProjectionAlgorithms (ptr): (%p)\n", self );
	_Mesh_Algorithms_Print( self, stream );
}

void _Mesh_ProjectionAlgorithms_AssignFromXML( void* algorithms, Stg_ComponentFactory* cf, void* data ) {
	_Mesh_Algorithms_AssignFromXML( algorithms, cf, data );
   _Mesh_ProjectionAlgorithms_Init( algorithms );
}

void _Mesh_ProjectionAlgorithms_Build( void* algorithms, void* data ) {
    _Mesh_Algorithms_Build( algorithms, data );
}

void _Mesh_ProjectionAlgorithms_Initialise( void* algorithms, void* data ) {
   Mesh_ProjectionAlgorithms* self = (Mesh_ProjectionAlgorithms*)algorithms;

   _Mesh_Algorithms_Initialise( algorithms, data );

   if( self->mesh->elTypes[0]->elementHasPointFunc == FeMesh_ElementType_ElementHasPoint ) {
      if( Stg_Class_CheckType(self->mesh->elTypes[0]->mesh->generator, ProjectionGenerator_Type ) ) {
	 /* only use ray cast if it uses the ProjectionGenerator, because that make linear elements. 
	    Ray cast only works for linear guys so far */
	 self->mesh->elTypes[0]->elementHasPointFunc = Mesh_ProjectionAlgorithms_ElementHasPoint;
      }
      else {
	 self->mesh->elTypes[0]->elementHasPointFunc = Mesh_ProjectionAlgorithms_ElementHasPointGeneral;
      }
   }
}

void _Mesh_ProjectionAlgorithms_Execute( void* algorithms, void* data ) {
    _Mesh_Algorithms_Execute( algorithms, data );
}

void _Mesh_ProjectionAlgorithms_Destroy( void* algorithms, void* data ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;

	Mesh_ProjectionAlgorithms_Destruct( self );

   _Mesh_Algorithms_Destroy( algorithms, data );
}

void Mesh_ProjectionAlgorithms_SetMesh( void* algorithms, void* mesh ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_ProjectionAlgorithms ) );

	Mesh_ProjectionAlgorithms_Destruct( self );
	_Mesh_Algorithms_SetMesh( self, mesh );
}

void Mesh_ProjectionAlgorithms_Update( void* algorithms ) {
	Mesh_ProjectionAlgorithms*	self = (Mesh_ProjectionAlgorithms*)algorithms;

	assert( self && Stg_CheckType( self, Mesh_ProjectionAlgorithms ) );
	assert( self->mesh );

	Mesh_ProjectionAlgorithms_Destruct( self );
	_Mesh_Algorithms_Update( self );

}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

Bool Mesh_ProjectionAlgorithms_ElementHasPointGeneral(
   void*             hexType,
   unsigned          elInd,
   double*           point, 
   MeshTopology_Dim* dim,
   unsigned*         ind )
{
   FeMesh              *mesh=NULL;
   const double        epsilon = 1e-7;
   FeMesh_ElementType* self = (FeMesh_ElementType*)hexType;
   int                 nDims, ii;
   Bool                border = False;

   IArray              *inc=NULL;
   ElementType         *elType=NULL;

   assert( self && Stg_CheckType( self, FeMesh_ElementType ) );

   mesh = (FeMesh*)self->mesh;

   // if not in cell, leave algorithm immediately
   if( !RayCast_IsPointInFace( mesh, elInd, point ) ) {
      return False;
   }

   elType = FeMesh_GetElementType( mesh, elInd );
   nDims = Mesh_GetDimSize( mesh );
   assert( nDims <= 3 );

   for( ii = 0; ii < nDims; ii++ ) {
      if( self->local[ii] < -1.0 - epsilon ||
         self->local[ii] > 1.0 + epsilon ) {
         return False;
      }

      if( self->local[ii] < -1.0 + epsilon ||
         self->local[ii] > 1.0 - epsilon ) {
         border = True;
      }
   }

   if( border ) {
      // Find the lowest indexed owning element.
      inc = IArray_New();
      double local[3];
      int    myGlobal, otherGlobal;
      int    jj, kk;

      myGlobal = FeMesh_ElementDomainToGlobal( mesh, elInd );

      Mesh_GetIncidence( mesh, nDims, elInd, nDims, inc );
      for( jj = 0; jj < inc->size; jj++ ) {
         FeMesh_CoordGlobalToLocal( mesh, inc->ptr[jj], point, local );
         for( kk = 0; kk < nDims; kk++ ) {
            if( local[kk] < -1.0 - epsilon ||
               local[kk] > 1.0 + epsilon ) {
               // This guy does not want the point.
               break;
            }
         }
         if( kk == nDims ) {
            // This guy wants the point too. Lowest index
            // wins.
            otherGlobal = FeMesh_ElementDomainToGlobal( mesh, inc->ptr[jj] );
            if( otherGlobal < myGlobal ) {
               Stg_Class_Delete( inc );
               return False;
            }
         }
      }
      Stg_Class_Delete( inc );
   }
	
   /* sanity test if node coord agree with local coords */
   int n_i;
   int *nodes, nNodes;
   double shapeFuncs[4], coord[3], *nodeCoord, lengthScale;

   inc = IArray_New();
   Mesh_GetIncidence( mesh, Mesh_GetDimSize(mesh), elInd, MT_VERTEX, inc );
   //Mesh_GetIncidence( mesh, nDims, elInd, 0, inc );
   nNodes = IArray_GetSize( inc );
   nodes = IArray_GetPtr( inc );

   memset( coord, 0, sizeof(double)*3);
   ElementType_EvaluateShapeFunctionsAt( elType, self->local, shapeFuncs );
   for( n_i=0; n_i<nNodes; n_i++ ) {
      nodeCoord = Mesh_GetVertex(mesh,nodes[n_i]); // get node coords

      // interpolate with shapefuncs node coord to local coord
      coord[0] = coord[0] + shapeFuncs[n_i]*nodeCoord[0];
      coord[1] = coord[1] + shapeFuncs[n_i]*nodeCoord[1];
      coord[2] = coord[2] + shapeFuncs[n_i]*nodeCoord[2];
   }

   // calc a length scale for error measure
   lengthScale = fabs( (Mesh_GetVertex(mesh, nodes[0])[0]) - (Mesh_GetVertex(mesh, nodes[3])[0]) ) ;
   assert( lengthScale > DBL_MIN );

   for( n_i = 0; n_i<nDims; n_i++ ) {
      if( fabs(coord[n_i] - point[n_i]) > 1e-7*lengthScale ) {
	 Stg_Class_Delete( inc );
	 return False;
      }	
   }

   Stg_Class_Delete( inc );
   *dim = nDims;
   *ind = elInd;
   return True;
}



Bool Mesh_ProjectionAlgorithms_ElementHasPoint(
   void*             hexType,
   unsigned          elInd,
   double*           point, 
   MeshTopology_Dim* dim,
   unsigned*         ind )
{
   FeMesh              *mesh=NULL;
   const double        epsilon = 1e-7;
   FeMesh_ElementType* self = (FeMesh_ElementType*)hexType;
   int                 nDims, ii;
   Bool                border = False;

   IArray              *inc=NULL;
   ElementType         *elType=NULL;

   assert( self && Stg_CheckType( self, FeMesh_ElementType ) );

   mesh = (FeMesh*)self->mesh;

   // if not in cell, leave algorithm immediately
   if( !RayCast_IsPointInFace( mesh, elInd, point ) ) {
      return False;
   }

   elType = FeMesh_GetElementType( mesh, elInd );
   nDims = Mesh_GetDimSize( mesh );
   assert( nDims <= 3 );

   // generate possible local coord
   FeMesh_CoordGlobalToLocal( mesh, elInd, point, self->local );

   for( ii = 0; ii < nDims; ii++ ) {
      if( self->local[ii] < -1.0 - epsilon ||
         self->local[ii] > 1.0 + epsilon ) {
         return False;
      }

      if( self->local[ii] < -1.0 + epsilon ||
         self->local[ii] > 1.0 - epsilon ) {
         border = True;
      }
   }

   if( border ) {
      // Find the lowest indexed owning element.
      inc = IArray_New();
      double local[3];
      int    myGlobal, otherGlobal;
      int    jj, kk;

      myGlobal = FeMesh_ElementDomainToGlobal( mesh, elInd );

      Mesh_GetIncidence( mesh, nDims, elInd, nDims, inc );
      for( jj = 0; jj < inc->size; jj++ ) {
         FeMesh_CoordGlobalToLocal( mesh, inc->ptr[jj], point, local );
         for( kk = 0; kk < nDims; kk++ ) {
            if( local[kk] < -1.0 - epsilon ||
               local[kk] > 1.0 + epsilon ) {
               // This guy does not want the point.
               break;
            }
         }
         if( kk == nDims ) {
            // This guy wants the point too. Lowest index
            // wins.
            otherGlobal = FeMesh_ElementDomainToGlobal( mesh, inc->ptr[jj] );
            if( otherGlobal < myGlobal ) {
               Stg_Class_Delete( inc );
               return False;
            }
         }
      }
      Stg_Class_Delete( inc );
   }
	
#if 0
   /* sanity test if node coord agree with local coords */
   int n_i;
   int *nodes, nNodes;
   double shapeFuncs[4], coord[3], *nodeCoord, lengthScale;

   inc = IArray_New();
   Mesh_GetIncidence( mesh, Mesh_GetDimSize(mesh), elInd, MT_VERTEX, inc );
   //Mesh_GetIncidence( mesh, nDims, elInd, 0, inc );
   nNodes = IArray_GetSize( inc );
   nodes = IArray_GetPtr( inc );

   memset( coord, 0, sizeof(double)*3);
   ElementType_EvaluateShapeFunctionsAt( elType, self->local, shapeFuncs );
   for( n_i=0; n_i<nNodes; n_i++ ) {
      nodeCoord = Mesh_GetVertex(mesh,nodes[n_i]); // get node coords

      // interpolate with shapefuncs node coord to local coord
      coord[0] = coord[0] + shapeFuncs[n_i]*nodeCoord[0];
      coord[1] = coord[1] + shapeFuncs[n_i]*nodeCoord[1];
      coord[2] = coord[2] + shapeFuncs[n_i]*nodeCoord[2];
   }

   // calc a length scale for error measure
   lengthScale = fabs( (Mesh_GetVertex(mesh, nodes[0])[0]) - (Mesh_GetVertex(mesh, nodes[3])[0]) ) ;
   assert( lengthScale > DBL_MIN );

   for( n_i = 0; n_i<nDims; n_i++ ) {
      if( fabs(coord[n_i] - point[n_i]) > 1e-7*lengthScale ) {
	 Stg_Class_Delete( inc );
	 return False;
      }	
   }

   Stg_Class_Delete( inc );
#endif
   *dim = nDims;
   *ind = elInd;
   return True;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Mesh_ProjectionAlgorithms_Destruct( Mesh_ProjectionAlgorithms* self ) {
	assert( self && Stg_CheckType( self, Mesh_ProjectionAlgorithms ) );
}


