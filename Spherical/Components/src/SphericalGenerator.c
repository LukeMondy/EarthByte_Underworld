#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"

#define OFFSET_TAG 6

/* Textual name of this class */
const Type SphericalGenerator_Type = "SphericalGenerator";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

void* _SphericalGenerator_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                                    _sizeOfSelf = sizeof(SphericalGenerator);
   Type                                                            type = SphericalGenerator_Type;
   Stg_Class_DeleteFunction*                                    _delete = _CartesianGenerator_Delete;
   Stg_Class_PrintFunction*                                      _print = _CartesianGenerator_Print;
   Stg_Class_CopyFunction*                                        _copy = NULL;
   Stg_Component_DefaultConstructorFunction*        _defaultConstructor = _SphericalGenerator_DefaultNew;
   Stg_Component_ConstructFunction*                          _construct = _SphericalGenerator_AssignFromXML;
   Stg_Component_BuildFunction*                                  _build = _CartesianGenerator_Build;
   Stg_Component_InitialiseFunction*                        _initialise = _SphericalGenerator_Initialise;
   Stg_Component_ExecuteFunction*                              _execute = _CartesianGenerator_Execute;
   Stg_Component_DestroyFunction*                              _destroy = _CartesianGenerator_Destroy;
   AllocationType                                    nameAllocationType = NON_GLOBAL;
   MeshGenerator_SetDimSizeFunc*                         setDimSizeFunc = CartesianGenerator_SetDimSize;
   MeshGenerator_GenerateFunc*                             generateFunc = CartesianGenerator_Generate;
   CartesianGenerator_SetTopologyParamsFunc*      setTopologyParamsFunc = _CartesianGenerator_SetTopologyParams;
   CartesianGenerator_GenElementsFunc*                  genElementsFunc = _CartesianGenerator_GenElements;
   CartesianGenerator_GenFacesFunc*                        genFacesFunc = _CartesianGenerator_GenFaces;
   CartesianGenerator_GenEdgesFunc*                        genEdgesFunc = _CartesianGenerator_GenEdges;
   CartesianGenerator_GenVerticesFunc*                  genVerticesFunc = _CartesianGenerator_GenVertices;
   CartesianGenerator_GenElementVertexIncFunc*  genElementVertexIncFunc = _CartesianGenerator_GenElementVertexInc;
   CartesianGenerator_GenVolumeEdgeIncFunc*        genVolumeEdgeIncFunc = _CartesianGenerator_GenVolumeEdgeInc;
   CartesianGenerator_GenVolumeFaceIncFunc*        genVolumeFaceIncFunc = _CartesianGenerator_GenVolumeFaceInc;
   CartesianGenerator_GenFaceVertexIncFunc*        genFaceVertexIncFunc = _CartesianGenerator_GenFaceVertexInc;
   CartesianGenerator_GenFaceEdgeIncFunc*          genFaceEdgeIncFunc = _CartesianGenerator_GenFaceEdgeInc;
   CartesianGenerator_GenEdgeVertexIncFunc*        genEdgeVertexIncFunc = _CartesianGenerator_GenEdgeVertexInc;
   CartesianGenerator_GenElementTypesFunc*         genElementTypesFunc = _SphericalGenerator_GenElementTypes;
   CartesianGenerator_CalcGeomFunc*                calcGeomFunc = SphericalGenerator_CalcCurvilinearGeom;

   SphericalGenerator* self = _SphericalGenerator_New(  SPHERICALGENERATOR_PASSARGS  );

   return self;
}

SphericalGenerator* _SphericalGenerator_New(  SPHERICALGENERATOR_DEFARGS  )
{
   SphericalGenerator* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(SphericalGenerator) );
   self = (SphericalGenerator*)_CartesianGenerator_New(  CARTESIANGENERATOR_PASSARGS  );

   /* Virtual info */
   self->setTopologyParamsFunc = setTopologyParamsFunc;
   self->genElementsFunc = genElementsFunc;
   self->genFacesFunc = genFacesFunc;
   self->genEdgesFunc = genEdgesFunc;
   self->genVerticesFunc = genVerticesFunc;
   self->genElementVertexIncFunc = genElementVertexIncFunc;
   self->genVolumeEdgeIncFunc = genVolumeEdgeIncFunc;
   self->genVolumeFaceIncFunc = genVolumeFaceIncFunc;
   self->genFaceVertexIncFunc = genFaceVertexIncFunc;
   self->genFaceEdgeIncFunc = genFaceEdgeIncFunc;
   self->genEdgeVertexIncFunc = genEdgeVertexIncFunc;
   self->genElementTypesFunc = genElementTypesFunc;
   self->calcGeomFunc = calcGeomFunc;


   return self;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SphericalGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data )
{
   SphericalGenerator*		self = (SphericalGenerator*)meshGenerator;

   _CartesianGenerator_AssignFromXML( meshGenerator, cf, data );

   // check domain size is valid
   Journal_Firewall( !(self->crdMin[0] <= 0 || self->crdMax[0] < self->crdMin[0]), global_error_stream,
        "Error in %s: Radius definition is wrong. Ensure maxX > minX & minX > 0\n", __func__ );

   // theta = [-180,180],
   Journal_Firewall(!(self->crdMin[1] < -180 || self->crdMax[1] > 180 || self->crdMax[1] < self->crdMin[1] ),
         global_error_stream,
         "Error in %s: Theta definition is wrong. Ensure -180 <= minY < maxY < 180\n", __func__ );

   if( self->nDims > 2 ) {
      // phi = (pi/2, pi)
      Journal_Firewall( !(self->crdMin[2] < -90 || self->crdMin[2] > 90 || self->crdMax[2] < self->crdMin[2]),
         global_error_stream, "\nError in %s: Phi definition is wrong. Ensure -90 < minZ < maxZ < 90\n", __func__);
   }

   self->fullAnnulus = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"FullAnnulus", False  );
   if( self->fullAnnulus )
   {
      self->crdMin[1] = -180;
      self->crdMax[1] = 180;
   }
}

void _SphericalGenerator_Initialise( void* meshGenerator, void* data )
{
   SphericalGenerator*		self = (SphericalGenerator*)meshGenerator;

   _CartesianGenerator_Initialise( meshGenerator, data );

}

void _SphericalGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh )
{
   SphericalGenerator*	self = (SphericalGenerator*)meshGenerator;
   Stream*			stream;
   unsigned		nDomainEls;
   unsigned		e_i;

   assert( self && Stg_CheckType( self, SphericalGenerator ) );

   stream = Journal_Register( Info_Type, (Name)self->type  );
   Journal_Printf( stream, "Generating element types...\n" );
   Stream_Indent( stream );

   mesh->nElTypes = 1;
   mesh->elTypes = Memory_Alloc_Array( Mesh_ElementType*, mesh->nElTypes, "Mesh::elTypes" );
   mesh->elTypes[0] = (Mesh_ElementType*)Mesh_HexType_New();
   Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
   nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
   mesh->elTypeMap = Memory_Alloc_Array( unsigned, nDomainEls, "Mesh::elTypeMap" );
   for( e_i = 0; e_i < nDomainEls; e_i++ )
      mesh->elTypeMap[e_i] = 0;

   Mesh_SetAlgorithms( mesh, Mesh_SphericalAlgorithms_New( "", NULL ) );

   MPI_Barrier( self->mpiComm );
   Journal_Printf( stream, "... element types are '%s',\n", mesh->elTypes[0]->type );
   Journal_Printf( stream, "... mesh algorithm type is '%s',\n", mesh->algorithms->type );
   Journal_Printf( stream, "... done.\n" );
   Stream_UnIndent( stream );
}

void SphericalGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
   /* function for creating a curvilinear mesh */
   SphericalGenerator* self = (SphericalGenerator*)_self;

   unsigned		n_i, d_i;
   double*         	vert;
   unsigned        	gNode;

   double dimRes[3];
   double radius,theta;

   /*
    * grid->sizes[] ... an array containing the number of nodes in a direction
    * gNode         ... global node id
    * inds          ... an ijk index to the node
    * steps         ... an array containing the physics length in each dim
    */
   for( d_i = 0 ; d_i < mesh->topo->nDims ; d_i++ )
   {
      dimRes[d_i] = steps[d_i]/(double)(grid->sizes[d_i]-1);
   }
   memcpy( self->sph_res, dimRes, 3*sizeof(double) );

   if( mesh->topo->nDims == 2 )
   {
      /* Loop over domain nodes. */
      for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ )
      {
         gNode = Sync_DomainToGlobal( sync, n_i );
         Grid_Lift( grid, gNode, inds );
         vert = Mesh_GetVertex( mesh, n_i );

         radius = self->crdMin[0] + dimRes[0] * (double)inds[0];
         theta = (M_PI/180) * (self->crdMin[1]+dimRes[1]*(double)inds[1]); // longitude position discretization

         vert[0] = radius * cos(theta);
         vert[1] = radius * sin(theta);
      }
   }
   else
   {
      double phi;
      /* Loop over domain nodes. */
      for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ )
      {
         vert = Mesh_GetVertex( mesh, n_i );

         gNode = Sync_DomainToGlobal( sync, n_i );
         Grid_Lift( grid, gNode, inds );

         radius = self->crdMin[0] + dimRes[0]*(double)inds[0];
         theta = (M_PI/180.0) * (self->crdMin[1]+dimRes[1]*(double)inds[1]); // longitude position discretization
         phi =   (M_PI/180.0) * (self->crdMin[2]+dimRes[2]*(double)inds[2]); // latitude position discretization

         vert[0] = radius * cos(theta) * cos(phi);
         vert[1] = radius * sin(phi);
         vert[2] = -radius * sin(theta) * cos(phi);
      }
   }
}
