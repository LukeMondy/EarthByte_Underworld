#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"
#include "Q2SphericalGenerator.h"


/* Textual name of this class */
const Type Q2SphericalGenerator_Type = "Q2SphericalGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Q2SphericalGenerator* Q2SphericalGenerator_New( Name name, AbstractContext* context )
{
   /* Variables set in this function */
   SizeT                                                    _sizeOfSelf = sizeof(Q2SphericalGenerator);
   Type                                                            type = Q2SphericalGenerator_Type;
   Stg_Class_DeleteFunction*                                    _delete = _C2Generator_Delete;
   Stg_Class_PrintFunction*                                      _print = _C2Generator_Print;
   Stg_Class_CopyFunction*                                        _copy = NULL;
   Stg_Component_DefaultConstructorFunction*        _defaultConstructor = (void* (*)(Name))_Q2SphericalGenerator_New;
   Stg_Component_ConstructFunction*                          _construct = _SphericalGenerator_AssignFromXML;
   Stg_Component_BuildFunction*                                  _build = _C2Generator_Build;
   Stg_Component_InitialiseFunction*                        _initialise = _C2Generator_Initialise;
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
   CartesianGenerator_GenElementTypesFunc*          genElementTypesFunc = Q2SphericalGenerator_GenElementTypes;
   CartesianGenerator_CalcGeomFunc*                calcGeomFunc = Q2SphericalGenerator_CalcCurvilinearGeom;

   Q2SphericalGenerator* self = _Q2SphericalGenerator_New(  C2SPHERICALGENERATOR_PASSARGS  );

   _MeshGenerator_Init( (MeshGenerator*)self, context );
   _CartesianGenerator_Init( (CartesianGenerator*)self );
   _C2Generator_Init( self );
   return self;
}

Q2SphericalGenerator* _Q2SphericalGenerator_New(  C2SPHERICALGENERATOR_DEFARGS  )
{
   Q2SphericalGenerator*	self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(Q2SphericalGenerator) );
   self = (Q2SphericalGenerator*)_CartesianGenerator_New(  CARTESIANGENERATOR_PASSARGS  );

   return self;
}

void _Q2SphericalGenerator_Delete( void* generator )
{
   _C2Generator_Delete( generator );
}

void Q2SphericalGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh )
{
   Q2SphericalGenerator*	self = (Q2SphericalGenerator*)meshGenerator;
   Stream*		stream;
   unsigned	nDomainEls;
   unsigned	vertMap[8] = {0, 2, 6, 8, 18, 20, 24, 26};
   unsigned	e_i;

   assert( self );

   stream = Journal_Register( Info_Type, (Name)self->type  );
   Journal_Printf( stream, "Generating element types...\n" );
   Stream_Indent( stream );

   mesh->nElTypes = 1;
   mesh->elTypes = AllocArray( Mesh_ElementType*, mesh->nElTypes );
   mesh->elTypes[0] = (Mesh_ElementType*)Mesh_HexType_New();
   Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
   Mesh_HexType_SetVertexMap( mesh->elTypes[0], vertMap );
   nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
   mesh->elTypeMap = AllocArray( unsigned, nDomainEls );
   for( e_i = 0; e_i < nDomainEls; e_i++ )
      mesh->elTypeMap[e_i] = 0;

   Mesh_SetAlgorithms( mesh, Mesh_SphericalAlgorithms_New( "", NULL ) );

   MPI_Barrier( self->mpiComm );
   Journal_Printf( stream, "... element types are '%s',\n", mesh->elTypes[0]->type );
   Journal_Printf( stream, "... mesh algorithm type is '%s',\n", mesh->algorithms->type );
   Journal_Printf( stream, "... done.\n" );
   Stream_UnIndent( stream );
}

void Q2SphericalGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
   /* function for creating a curvilinear mesh */
   SphericalGenerator* self = (SphericalGenerator*)_self;

   unsigned		n_i, d_i;
   double*         	vert;
   unsigned        	gNode;

   double nodalRes[3];
   double elRes[3];
   double radius,theta;

   /*
    * grid->sizes[] ... an array containing the number of nodes in a direction
    * gNode         ... global node id
    * inds          ... an ijk index to the node
    * steps         ... an array containing the physics length in each dim
    */
   for( d_i = 0 ; d_i < mesh->topo->nDims ; d_i++ )
   {
      nodalRes[d_i] = steps[d_i]/(double)(grid->sizes[d_i]-1);
   }

   for( d_i = 0 ; d_i < mesh->topo->nDims ; d_i++ )
   {
      elRes[d_i] = steps[d_i]/(0.5*(double)(grid->sizes[d_i]-1));
   }
   memcpy( self->sph_res, elRes, 3*sizeof(double) );

   printf(" my dim res is (%g, %g, %g)\n", nodalRes[0], nodalRes[1], nodalRes[2]);
   /* Loop over domain nodes. */
   for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ )
   {
      gNode = Sync_DomainToGlobal( sync, n_i );
      Grid_Lift( grid, gNode, inds );
      vert = Mesh_GetVertex( mesh, n_i );

      radius = self->crdMin[0] + nodalRes[0] * (double)inds[0];
      theta = (M_PI/180) * (self->crdMin[1] + nodalRes[1]*(double)inds[1]); // longitude dis

      vert[0] = radius * cos(theta);
      vert[1] = radius * sin(theta);


   }
}
