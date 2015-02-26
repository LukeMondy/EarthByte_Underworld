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
const Type Q2RSGenerator_Type = "Q2RSGenerator";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

void* _Q2RSGenerator_DefaultNew( Name name, AbstractContext* context )
{
   /* Variables set in this function */
   SizeT                                                    _sizeOfSelf = sizeof(Q2RSGenerator);
   Type                                                            type = Q2RSGenerator_Type;
   Stg_Class_DeleteFunction*                                    _delete = _C2Generator_Delete;
   Stg_Class_PrintFunction*                                      _print = _C2Generator_Print;
   Stg_Class_CopyFunction*                                        _copy = NULL;
   Stg_Component_DefaultConstructorFunction*        _defaultConstructor = _Q2RSGenerator_DefaultNew;
   Stg_Component_ConstructFunction*                          _construct = _Q2RSGenerator_AssignFromXML;
   Stg_Component_BuildFunction*                                  _build = _C2Generator_Build;
   Stg_Component_InitialiseFunction*                        _initialise = _Q2RSGenerator_Initialise;
   Stg_Component_ExecuteFunction*                              _execute = _C2Generator_Execute;
   Stg_Component_DestroyFunction*                              _destroy = _C2Generator_Destroy;
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
   CartesianGenerator_GenElementTypesFunc*          genElementTypesFunc = _Q2RSGenerator_GenElementTypes;
   CartesianGenerator_CalcGeomFunc*                        calcGeomFunc = Q2RSGenerator_CalcCurvilinearGeom;

   Q2RSGenerator* self = _Q2RSGenerator_New( Q2RSGENERATOR_PASSARGS );

   _MeshGenerator_Init( (MeshGenerator*)self, context );
   _CartesianGenerator_Init( (CartesianGenerator*)self );
   _C2Generator_Init( self );

   return self;
}

Q2RSGenerator* _Q2RSGenerator_New( Q2RSGENERATOR_DEFARGS )
{
   Q2RSGenerator* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(Q2RSGenerator) );
   self = (Q2RSGenerator*)_CartesianGenerator_New( RSGENERATOR_PASSARGS );

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

void _Q2RSGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data )
{
   Q2RSGenerator*			self 		= (Q2RSGenerator*)meshGenerator;
   FiniteElementContext*	context 	= Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, NULL  );
   Dictionary*			dict		= Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, (Dictionary_Entry_Key)self->name ) );
   Dictionary_Entry_Value	*minList, *maxList, *tmp;
   char*			rootKey;
   unsigned			d_i;
   double 			maxVal;

   _CartesianGenerator_AssignFromXML( meshGenerator, cf, data );

   // mesh min/max coords are written to file from the mesh->min/maxGlobalCoord field. these differ from the min/max
   // coords as defined in the Q2RSGenerator->min/maxCrd (these are in r-theta-phi, not x-y-z), so if we're restarting
   // we need to over-ride the coords from the hdf5 mesh file with those from the xml (assuming they don't differ
   // between runs...
   if( self->nDims > 2 ) {
      if( context->loadFromCheckPoint ) {
         minList = Dictionary_Get( dict, (Dictionary_Entry_Key)"minCoord" );
         maxList = Dictionary_Get( dict, (Dictionary_Entry_Key)"maxCoord" );
         if( minList && maxList ) {
                assert( Dictionary_Entry_Value_GetCount( minList ) >= self->nDims );
                assert( Dictionary_Entry_Value_GetCount( maxList ) >= self->nDims  );
                for( d_i = 0; d_i < self->nDims; d_i++ ) {
                   tmp = Dictionary_Entry_Value_GetElement( minList, d_i );
                   rootKey = Dictionary_Entry_Value_AsString( tmp );

                   if( !Stg_StringIsNumeric( (char*)rootKey ) )
                      tmp = Dictionary_Get( cf->rootDict, (Dictionary_Entry_Key)rootKey );
                   self->crdMin[d_i] = Dictionary_Entry_Value_AsDouble( tmp );

                   tmp = Dictionary_Entry_Value_GetElement( maxList, d_i );
                   rootKey = Dictionary_Entry_Value_AsString( tmp );

                   if( !Stg_StringIsNumeric( (char*)rootKey ) )
                      tmp = Dictionary_Get( cf->rootDict, (Dictionary_Entry_Key)rootKey );
                   self->crdMax[d_i] = Dictionary_Entry_Value_AsDouble( tmp );
                   /* test to ensure provided domain is valid */
                   maxVal = (abs(self->crdMax[d_i]) > abs(self->crdMin[d_i])) ? abs(self->crdMax[d_i]) : abs(self->crdMin[d_i]);
                   if( maxVal == 0 ) maxVal = 1;  /* if maxVal is zero, then both numbers must be zero, set to one as next test will fail */
                   Journal_Firewall( ( ( (self->crdMax[d_i] - self->crdMin[d_i])/maxVal) > 1E-10 || d_i==J_AXIS), global_error_stream,
                     "\n\nError in %s for %s '%s'\n\n"
                     "Dimension of domain (min = %f, max = %f) for component number %u is not valid.\n\n",
                     __func__, self->type, self->name, self->crdMin[d_i], self->crdMax[d_i], d_i );
                }
         }
      }
      // phi = (pi/2, pi)
      Journal_Firewall( !(fabs(self->crdMin[2] -self->crdMax[2])>=179.99 || self->crdMax[2] < self->crdMin[2]),
         global_error_stream, "\nError in %s: Phi definition is wrong. Ensure minZ < maxZ & abs(minZ-maxZ)<180\n", __func__);
   }
   Journal_Firewall( self->nDims==3, global_error_stream,
        "Error in %s: Must have 3 dimensions for component %s\n", __func__, self->name );
   // check domain size is valid
   Journal_Firewall( !(self->crdMin[0] <= 0 || self->crdMax[0] < self->crdMin[0]), global_error_stream,
        "Error in %s: Radius definition is wrong. Ensure maxX > minX & minX > 0\n", __func__ );

   Journal_Firewall(!( fabs(self->crdMin[1]-self->crdMax[1]) >= 179.99 || self->crdMax[1] < self->crdMin[1] ),
         global_error_stream, "Error in %s: Theta definition is wrong. Ensure minY < maxY && abs(minY-maxY) < 180\n", __func__ );
}

void _Q2RSGenerator_Initialise( void* meshGenerator, void* data )
{
   Q2RSGenerator*		self = (Q2RSGenerator*)meshGenerator;

   _CartesianGenerator_Initialise( meshGenerator, data );
   _C2Generator_Initialise( meshGenerator, data );
}

void _Q2RSGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh )
{
   Q2RSGenerator*	self = (Q2RSGenerator*)meshGenerator;
   Stream*			stream;
   unsigned		nDomainEls;
   unsigned		e_i;

   assert( self && Stg_CheckType( self, Q2RSGenerator ) );

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

void Q2RSGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
   /* function for creating a curvilinear mesh */
   Q2RSGenerator* self = (Q2RSGenerator*)_self;

   unsigned		n_i, d_i;
   double*         	vert;
   unsigned        	gNode;

   double dimRes[3];
   double theta;

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

   double phi;
   double X,Y,r,d;

   /* Loop over domain nodes. */
   for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ )
   {
      vert = Mesh_GetVertex( mesh, n_i );

      gNode = Sync_DomainToGlobal( sync, n_i );
      Grid_Lift( grid, gNode, inds );

      r = self->crdMin[0] + dimRes[0]*(double)inds[0];
      theta = (M_PI/180.0) * (self->crdMin[1]+dimRes[1]*(double)inds[1]); // longitude position discretization
      phi =   (M_PI/180.0) * (self->crdMin[2]+dimRes[2]*(double)inds[2]); // latitude position discretization

      X = tan(theta);
      Y = tan(phi);
      d = sqrt(1 + X*X + Y*Y);

      vert[0] = r/d * X;
      vert[1] = r/d * Y;
      vert[2] = r/d * 1;
   }
}
