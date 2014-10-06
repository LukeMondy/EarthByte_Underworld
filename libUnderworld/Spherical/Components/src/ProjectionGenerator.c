#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"

#define OFFSET_TAG 6

/* Textual name of this class */
const Type ProjectionGenerator_Type = "ProjectionGenerator";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

ProjectionGenerator* _ProjectionGenerator_DefaultNew( Name name )
{
  /* Variables set in this function */
  SizeT                                                    _sizeOfSelf = sizeof(ProjectionGenerator);
  Type                                                            type = ProjectionGenerator_Type;
  Stg_Class_DeleteFunction*                                    _delete = _CartesianGenerator_Delete;
  Stg_Class_PrintFunction*                                      _print = _CartesianGenerator_Print;
  Stg_Class_CopyFunction*                                        _copy = NULL;
  Stg_Component_DefaultConstructorFunction*        _defaultConstructor = (void* (*)(Name))_CartesianGenerator_New;
  Stg_Component_ConstructFunction*                          _construct = _ProjectionGenerator_AssignFromXML;
  Stg_Component_BuildFunction*                                  _build = _ProjectionGenerator_Build;
  Stg_Component_InitialiseFunction*                        _initialise = _ProjectionGenerator_Initialise;
  Stg_Component_ExecuteFunction*                              _execute = _CartesianGenerator_Execute;
  Stg_Component_DestroyFunction*                              _destroy = _ProjectionGenerator_Destroy;
  AllocationType                                    nameAllocationType = NON_GLOBAL;
  MeshGenerator_SetDimSizeFunc*                         setDimSizeFunc = NULL;//CartesianGenerator_SetDimSize;
  MeshGenerator_GenerateFunc*                             generateFunc = ProjectionGenerator_Generate;
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
  CartesianGenerator_GenElementTypesFunc*         genElementTypesFunc = _ProjectionGenerator_GenElementTypes;
  CartesianGenerator_CalcGeomFunc*                calcGeomFunc = NULL;

  ProjectionGenerator* self = _ProjectionGenerator_New(  PROJECTIONGENERATOR_PASSARGS  );

  return self;
}

void _ProjectionGenerator_Destroy( void* meshGenerator, void* data ) {
   ProjectionGenerator *self = (ProjectionGenerator*)meshGenerator;

}
ProjectionGenerator* _ProjectionGenerator_New(  PROJECTIONGENERATOR_DEFARGS  )
{
  ProjectionGenerator* self;

  /* Allocate memory */
  assert( _sizeOfSelf >= sizeof(ProjectionGenerator) );
  self = (ProjectionGenerator*)_CartesianGenerator_New(  CARTESIANGENERATOR_PASSARGS  );

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

void _ProjectionGenerator_Init( ProjectionGenerator* self, int nSquares, int map, double distance, double decay, double superel )
{
  unsigned size[3];

  assert( nSquares > 0 );
  assert( distance > 0 );
  assert( decay > 0 );

  //self->nSquares = nSquares;
  self->superel = superel;

  // reset the elementGrid data struct to contain the correct number of elements
  //size[0] = size[1] = size[2] = 2*nSquares;
  //Grid_SetSizes( self->elGrid, size );

  if( map == 1 )
  {
    self->_mapRing = ProjectionGenerator_ExpMap;
    assert( decay > 0 );
  }
  else
  {
    self->_mapRing = ProjectionGenerator_PowerMap;
  }

  self->distance = distance;
  self->decay = decay;

  if( self->nDims == 2 ) {
     if( !self->equiangle )
       self->calcGeomFunc = ProjectionGenerator_EquiDistanceGeom; // for equidistanc
     else
       self->calcGeomFunc = ProjectionGenerator_EquiAngleGeom; // for equiangle
  } else {
     self->calcGeomFunc = ProjectionGenerator_EquiAngle3DGeom;
  }
}

void _ProjectionGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data )
{
  ProjectionGenerator*		self = (ProjectionGenerator*)meshGenerator;
  Mesh*  mesh;
  int nSquares;
  double distance, decay, superel;
  int mapping;

  // set the communicator
   MeshGenerator_SetMPIComm( self, MPI_COMM_WORLD );

   //get the context
   self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", AbstractContext, False, data );
	if( !self->context  )
		self->context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  );

   // get the mesh
   mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"mesh", Mesh, False, data );
	if( mesh  )
		MeshGenerator_AddMesh( self, mesh );

	_MeshGenerator_SetDimSize( self, 3 );

   // enable topographic mappings
 //  MeshGenerator_SetIncidenceState( self, 0, 0, True ); // node-to-node
 //  MeshGenerator_SetIncidenceState( self, 0, 3, True ); // node-to-element
 //  MeshGenerator_SetIncidenceState( self, 3, 0, True ); // element-to-node
 //  MeshGenerator_SetIncidenceState( self, 3, 3, True ); // element-to-element

   self->comm = Comm_New();
}

void ProjectionGenerator_Generate( void* meshGenerator, void* data ) {
   
  ProjectionGenerator*		self = (ProjectionGenerator*)meshGenerator;
  IGraph *topo = self->meshes[0]->topo; // TODO I hate this backwards pointer crap

   int v_i, n_i, e_i, nVerts, nEls;
   FILE *file=NULL;
   int **e_n, **nbrs;
   double *verts;
   int *locals, bc;
   SphericalFeMesh* mesh = (SphericalFeMesh*)self->meshes[0];

   MeshTopology_SetComm( topo, self->comm );
   MeshTopology_SetNumDims( topo, 3 );

   /* read in node coords */
   file = fopen( "coords.dat", "r" );
   assert(file);
   fscanf( file, "%d\n", &nVerts ); 

   mesh->bndNodeSet = IndexSet_New( nVerts );

   mesh->innerSet = IndexSet_New( nVerts );
   mesh->outerSet = IndexSet_New( nVerts );

   // allocate node coord array and fill in
   verts = self->meshes[0]->vertices = Memory_Alloc_Array_Unnamed( double, 3*nVerts );
   for( v_i=0; v_i<nVerts; v_i++ ) {

      fscanf( file, "%lg %lg %lg %d\n", &(verts[v_i]), &(verts[v_i+1]), &(verts[v_i+2]), &bc );

      locals[v_i]=v_i;

      // if bc add it to the correct set
      if( bc == 1 ) {
         IndexSet_Add( mesh->innerSet, v_i );
      } else if (bc == 2) {
         IndexSet_Add( mesh->outerSet, v_i );
      }
   }
   fclose( file );

   // build set counts
   IndexSet_UpdateMembersCount( mesh->innerSet );
   IndexSet_UpdateMembersCount( mesh->outerSet );

	IGraph_SetElements( topo, 0, nVerts, locals );
   FreeArray(locals);

   /* read in element node mapping */
   file = fopen( "elements.dat", "r" );
   assert(file);
   fscanf( file, "%d\n", &nEls ); 

   // allocate node coord array and fill in
   e_n = Memory_Alloc_Array_Unnamed( int*, nEls );
	locals = Memory_Alloc_Array_Unnamed( unsigned, nEls );
   for( e_i=0; e_i<nEls; e_i++ ) {
      locals[e_i]=e_i;
      e_n[e_i] = Memory_Alloc_Array_Unnamed( int, 8 );

      /* below we reorder the node indices because in elements.dat the ording is VTK_HEXAHEDRA
         which is different to StgFEM. The follow ordering creates StgFEM compliant orderings */ 
      fscanf( file, "%d ", &(e_n[e_i][0]) );
      fscanf( file, "%d ", &(e_n[e_i][1]) );
      fscanf( file, "%d ", &(e_n[e_i][3]) );
      fscanf( file, "%d ", &(e_n[e_i][2]) );
      fscanf( file, "%d ", &(e_n[e_i][4]) );
      fscanf( file, "%d ", &(e_n[e_i][5]) );
      fscanf( file, "%d ", &(e_n[e_i][7]) );
      fscanf( file, "%d\n", &(e_n[e_i][6]) );
      /*
      for( n_i=0; n_i<7; n_i++ ) {
         fscanf( file, "%d ", &(e_n[e_i][n_i]) );
      }
      fscanf( file, "%d\n", &(e_n[e_i][n_i]) );
      */
   }
   fclose( file );

	IGraph_SetElements( topo, 3, nEls, locals );
   FreeArray(locals);

   /* Now everything is read in, we set up the mesh */
   for(e_i=0; e_i<nEls; e_i++ ) {
		IGraph_SetIncidence( topo, MT_VOLUME, e_i, MT_VERTEX, 8, e_n[e_i] );
   }

   // build the node-to-element mapping
   IGraph_InvertIncidence( topo, MT_VERTEX, MT_VOLUME );
   /*
   for( v_i=0; v_i<nVerts; v_i++ ) {
		IGraph_SetIncidence( topo, MT_VERTEX, v_i, MT_VERTEX, 6, nbr[v_i] );
   }
   */
  _ProjectionGenerator_GenElementTypes( self, self->meshes[0] );
}
void _ProjectionGenerator_Build( void* meshGenerator, void* data )
{
  ProjectionGenerator*		self = (ProjectionGenerator*)meshGenerator;

  IGraph *topo = self->meshes[0]->topo; // TODO I hate this backwards pointer crap

   MeshTopology_SetComm( topo, self->comm );
   MeshTopology_SetNumDims( topo, 3 );


}
void _ProjectionGenerator_Initialise( void* meshGenerator, void* data )
{
  //ProjectionGenerator*		self = (ProjectionGenerator*)meshGenerator;

  _CartesianGenerator_Initialise( meshGenerator, data );
}

void _ProjectionGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh )
{
  ProjectionGenerator*	self = (ProjectionGenerator*)meshGenerator;
  Stream*			stream;
  unsigned		nDomainEls;
  unsigned		e_i;

  assert( self && Stg_CheckType( self, ProjectionGenerator ) );

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

  Mesh_SetAlgorithms( mesh, Mesh_ProjectionAlgorithms_New( "", NULL ) );

  MPI_Barrier( self->mpiComm );
  Journal_Printf( stream, "... element types are '%s',\n", mesh->elTypes[0]->type );
  Journal_Printf( stream, "... mesh algorithm type is '%s',\n", mesh->algorithms->type );
  Journal_Printf( stream, "... done.\n" );
  Stream_UnIndent( stream );
}

double ProjectionGenerator_PowerMap( void* _self, int sq_i )
{
  ProjectionGenerator* self = (ProjectionGenerator*)_self;
  double x = sq_i/(double)self->nSquares;

  return self->distance * pow(x, (1.0/self->decay) );
}

double ProjectionGenerator_ExpMap( void* _self, int sq_i )
{
  ProjectionGenerator* self = (ProjectionGenerator*)_self;
  double x = sq_i/(double)self->nSquares; // x = [0,1]

  // self->distance scales the mapping
  // (1-exp( -self->decay * x)) / (1-exp(-self->decay)) maps x to x', where x' = [0,1]
  return self->distance  * (1-exp( -self->decay * x)) / (1-exp(-self->decay));
}

void ProjectionGenerator_EquiDistanceGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
  /* Purpose to creat a cube-sphere mesh like mesh */

  ProjectionGenerator* self = (ProjectionGenerator*)_self;

  double dimRes[3], centre[3], vec[3], *vert;

  unsigned gNode, n_i, n_sideNodes, nDomain;
  int ii,jj;
  int nSquares = self->nSquares;

  // get the local number of nodes
  nDomain = Mesh_GetDomainSize( mesh, MT_VERTEX );

  // calc the number of sides nodes on the largest square
  n_sideNodes = self->vertRange[0];

  // calculate the cube resolution for the initial cube placement & centre of cube
  for( ii = 0 ; ii < mesh->topo->nDims ; ii++ )
  {
    dimRes[ii] = steps[ii]/(double)(n_sideNodes-1);
    centre[ii] = self->crdMin[ii] + (double)steps[ii]/2.0;
  }
  memcpy( self->pgen_res, dimRes, mesh->topo->nDims*sizeof(double) );

  /* Put nodes into a regular cube. x=[0,1], y=[0,1], z=[0,1] */
  for( n_i = 0; n_i < nDomain ; n_i++ )
  {
    vert = Mesh_GetVertex( mesh, n_i );

    // get ijk parameterisation of node
    gNode = Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i );
    Grid_Lift( grid, gNode, inds );

    vert[0] = self->crdMin[0] + (double)inds[0]*dimRes[0];
    vert[1] = self->crdMin[1] + (double)inds[1]*dimRes[1];
  }

  double mag, factor;
  double distance, n, frac,  xyz[2], theta;
  double che = self->superel;

  /* assign each node's sq_id */
  int sq_i, end;
  unsigned ind[3];
  ind[2]=0;
  for( sq_i=nSquares; sq_i>=0 ; sq_i-- )
  {
    /* In ijk parameters -
    	 start from corner (sq_i, sq_i)
    	 traverse along j and tag node's sq_id
    	 then traverse along i and tag node's sq_id
     */

    end = nSquares - sq_i; // calc end index for this square

    ii = end;
    for( jj=end; jj<n_sideNodes-end; jj++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
    }

    ii = n_sideNodes-end-1;
    for( jj=end; jj<n_sideNodes-end; jj++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
    }

    jj=end;
    for( ii=end; ii<n_sideNodes-end; ii++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
    }
    jj=n_sideNodes-end-1;
    for( ii=end; ii<n_sideNodes-end; ii++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
    }
  }

  // projection time

  for( n_i=0; n_i<nDomain; n_i++ )
  {
    /* project each point along vec, where the distance is a function of sq_i */
    vert = Mesh_GetVertex( mesh, n_i );

    // calculate the distance between n1 and n2
    StGermain_VectorSubtraction( vec, vert, centre, 2 ) ;

    mag = sqrt( vec[0]*vec[0] + vec[1]*vec[1] );

    if( mag < 1e-10 && self->nodeSqId[n_i] == 0 )
    {
      // this is the center node, leave it in the center
      vert[0] = centre[0];
      vert[1] = centre[1];
      continue;
    }

    xyz[0] = fabs(vec[0]);
    xyz[1] = fabs(vec[1]);

    n = che - ((che-2)*self->nodeSqId[n_i]/nSquares);

    distance = self->_mapRing( self, self->nodeSqId[n_i] );
    frac = pow(distance,n) / ( pow( xyz[0],n) + pow( xyz[1],n) );
    factor = pow( frac, 1/n);


    if(!finite( factor ) )
    {
      assert(0);
    }


    // calculate projection factor
    //factor = self->_mapRing( self, self->nodeSqId[n_i] ) / mag;

    // project that node!
    vert[0] = centre[0] + factor*vec[0];
    vert[1] = centre[1] + factor*vec[1];

    if( self->rotateMesh )
    {
      /* rotate nodes !!! */

      theta = self->nodeSqId[n_i] * 5 * M_PI / 180;

      xyz[0] = cos(theta)*vert[0] + sin(theta)*vert[1];
      xyz[1] = -sin(theta)*vert[0] + cos(theta)*vert[1];

      vert[0] = xyz[0];
      vert[1] = xyz[1];
    }
  }
#if 0
#endif

}


void ProjectionGenerator_EquiAngle3DGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
  /* Purpose to creat a cube-sphere mesh like mesh */

   ProjectionGenerator* self = (ProjectionGenerator*)_self;

   unsigned		n_i, d_i;
   double*         	vert;
   unsigned        	gNode;

   double dimRes[3];
   double radius,theta;
   double phi=0;
   double X,Y,r,R;

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
   memcpy( self->pgen_res, dimRes, mesh->topo->nDims*sizeof(double) );


   /* Loop over domain nodes. */
   for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ )
   {
      gNode = Sync_DomainToGlobal( sync, n_i );
      Grid_Lift( grid, gNode, inds );
      vert = Mesh_GetVertex( mesh, n_i );

      theta = self->crdMin[0] + dimRes[0] * (double)inds[0];
      radius = self->crdMin[1] + dimRes[1] *(double)inds[1];
      phi = self->crdMin[2] + dimRes[2] *(double)inds[2];

      X = radius*tan(theta);
      Y = radius*tan(phi);
      r = sqrt(radius*radius + X*X + Y*Y);
      R = sqrt(3)*radius;

      vert[0] = R/r * X;
      vert[1] = R/r * Y;
      vert[2] = R/r * radius;

   }
}

void ProjectionGenerator_EquiAngleGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps )
{
  /* Purpose to creat a cube-sphere mesh like mesh */

  ProjectionGenerator* self = (ProjectionGenerator*)_self;

  double dimRes[3], centre[3], vec[3], *vert;

  unsigned gNode, n_i, n_sideNodes, nDomain;
  int ii,jj;
  int nSquares = self->nSquares;

  //create storage and calculate the delta theta for each ring of nodes
  double *dt = Memory_Alloc_Array_Unnamed( double, nSquares );
  for( ii=0; ii<nSquares; ii++ )
  {
    // 2*(ii+i)+1 is the number of nodes per side
    dt[ii] = 0.5*M_PI/(2*(ii+1)+1-1);
  }

  // get the local number of nodes
  nDomain = Mesh_GetDomainSize( mesh, MT_VERTEX );

  // calc the number of sides nodes on the largest square
  n_sideNodes = self->vertRange[0];

  // calculate the cube resolution for the initial cube placement & centre of cube
  for( ii = 0 ; ii < mesh->topo->nDims ; ii++ )
  {
    dimRes[ii] = 0.5*M_PI/(double)(n_sideNodes-1);
    centre[ii] = self->crdMin[ii] + (double)steps[ii]/2.0;
  }
   memcpy( self->pgen_res, dimRes, mesh->topo->nDims*sizeof(double) );

  double mag, factor;
  double distance, n, frac,  xyz[2], theta;
  double che = self->superel;

  /* assign each node's sq_id */
  int sq_i, end;
  unsigned ind[3];
  ind[2]=0;
  for( sq_i=nSquares; sq_i>=0 ; sq_i-- )
  {
    /* In ijk parameters -
    	 start from corner (sq_i, sq_i)
    	 traverse along j and tag node's sq_id
    	 then traverse along i and tag node's sq_id
     */

    end = nSquares - sq_i; // calc end index for this square
    distance = self->_mapRing( self, sq_i );
    theta = 5 *M_PI/4;

    // left wall
    ii = end;
    for( jj=end; jj<n_sideNodes-end; jj++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;

      vert = Mesh_GetVertex( mesh, n_i );
      vert[0] = distance * cos( theta );
      vert[1] = distance * sin( theta );
      printf("vert %d, ring %d, theta %g (%g, %g)\n", n_i, sq_i, theta, vert[0], vert[1] );
      theta = theta - dt[sq_i-1];
    }

    theta = theta + dt[sq_i-1]; //dodgy fix
    // top wall
    jj=n_sideNodes-end-1;
    for( ii=end; ii<n_sideNodes-end; ii++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
      vert = Mesh_GetVertex( mesh, n_i );
      vert[0] = distance * cos( theta );
      vert[1] = distance * sin( theta );
      printf("vert %d, ring %d, theta %g (%g, %g)\n", n_i, sq_i, theta, vert[0], vert[1] );
      theta = theta - dt[sq_i-1];
    }
    theta = theta + dt[sq_i-1]; //dodgy fix
    // right wall
    ii = n_sideNodes-end-1;
    //for( jj=end; jj<n_sideNodes-end; jj++ )
    for( jj=n_sideNodes-end-1; jj>=end; jj-- )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
      vert = Mesh_GetVertex( mesh, n_i );
      vert[0] = distance * cos( theta );
      vert[1] = distance * sin( theta );
      printf("vert %d, ring %d, theta %g (%g, %g)\n", n_i, sq_i, theta, vert[0], vert[1] );
      theta = theta - dt[sq_i-1];
    }

    theta = theta + dt[sq_i-1]; //dodgy fix
    // bottom wall
    jj=end;
    for( ii=n_sideNodes-end-1; ii>=end; ii-- )
      //for( ii=end; ii<n_sideNodes-end; ii++ )
    {
      ind[0]=ii;
      ind[1]=jj;
      Mesh_GlobalToDomain( mesh, MT_VERTEX, Grid_Project( grid, ind ), &n_i );
      // if node is local we set its square index
      if( n_i < nDomain ) self->nodeSqId[n_i] = sq_i;
      vert = Mesh_GetVertex( mesh, n_i );
      vert[0] = distance * cos( theta );
      vert[1] = distance * sin( theta );
      printf("vert %d, ring %d, theta %g (%g, %g)\n", n_i, sq_i, theta, vert[0], vert[1] );
      theta = theta - dt[sq_i-1];
    }
    theta = theta + dt[sq_i-1]; //dodgy fix

  }

  if( self->rotateMesh )
  {
    /* rotate nodes */
    for( n_i = 0; n_i < nDomain ; n_i++ )
    {
      vert = Mesh_GetVertex( mesh, n_i );

      /* rotate nodes !!! */
      theta =  M_PI / 4;

      xyz[0] = cos(theta)*vert[0] + sin(theta)*vert[1];
      xyz[1] = -sin(theta)*vert[0] + cos(theta)*vert[1];

      vert[0] = xyz[0];
      vert[1] = xyz[1];

    }
  }


  // projection time

#if 0
  double zero[2] = {0.0,0.0};
  for( n_i=0; n_i<nDomain; n_i++ )
  {
    /* project each point along vec, where the distance is a function of sq_i */
    vert = Mesh_GetVertex( mesh, n_i );

    // calculate the distance between n1 and n2
    StGermain_VectorSubtraction( vec, vert, zero, 2 ) ;

    mag = sqrt( vec[0]*vec[0] + vec[1]*vec[1] );

    if( mag < 1e-10 && self->nodeSqId[n_i] == 0 )
    {
      // this is the center node, leave it in the center
      vert[0] = centre[0];
      vert[1] = centre[1];
      continue;
    }

#if 0
    xyz[0] = fabs(vec[0]);
    xyz[1] = fabs(vec[1]);

    n = che - ((che-2)*self->nodeSqId[n_i]/nSquares);

    distance = self->_mapRing( self, self->nodeSqId[n_i] );
    frac = pow(distance,n) / ( pow( xyz[0],n) + pow( xyz[1],n) );
    factor = pow( frac, 1/n);
#endif

    // calculate projection factor
    factor = self->_mapRing( self, self->nodeSqId[n_i] ) / mag;

    if(!isfinite( factor ) )
    {
      assert(0);
    }

    // project that node!
    vert[0] = centre[0] + factor*vec[0];
    vert[1] = centre[1] + factor*vec[1];

  }
#endif

}
