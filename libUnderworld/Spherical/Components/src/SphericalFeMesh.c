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
const Type SphericalFeMesh_Type = "SphericalFeMesh";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

SphericalFeMesh* SphericalFeMesh_New( Name name, AbstractContext* context )
{
    SphericalFeMesh* self = _SphericalFeMesh_DefaultNew( name );
    _Mesh_Init( (Mesh*)self, context );
    /* SphericalFeMesh info */
    _SphericalFeMesh_Init( self, NULL, NULL, False ); /* this is a useless Init() */

    return self;
}

SphericalFeMesh* _SphericalFeMesh_DefaultNew( Name name )
{
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SphericalFeMesh);
    Type                                                      type = SphericalFeMesh_Type;
    Stg_Class_DeleteFunction*                              _delete = _SphericalFeMesh_Delete;
    Stg_Class_PrintFunction*                                _print = _SphericalFeMesh_Print;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (void* (*)(Name))_SphericalFeMesh_New;
    Stg_Component_ConstructFunction*                    _construct = _SphericalFeMesh_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SphericalFeMesh_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SphericalFeMesh_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SphericalFeMesh_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SphericalFeMesh_Destroy;
    AllocationType                              nameAllocationType = NON_GLOBAL;

    /* The following terms are parameters that have been passed into or defined in this function but are being set before being passed onto the parent */
    Stg_Class_CopyFunction*        _copy = NULL;

    return _SphericalFeMesh_New(  FEMESH_PASSARGS  );
}

SphericalFeMesh* _SphericalFeMesh_New(  FEMESH_DEFARGS  )
{
    SphericalFeMesh*	self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SphericalFeMesh) );
    self = (SphericalFeMesh*)_Mesh_New(  MESH_PASSARGS  );

    return self;
}

void _SphericalFeMesh_Init( SphericalFeMesh* self, ElementType* elType, const char* family, Bool elementMesh )
{
    Stream*	stream;

    assert( self && Stg_CheckType( self, SphericalFeMesh ) );

    stream = Journal_Register( Info_Type, (Name)self->type  );
    Stream_SetPrintingRank( stream, 0 );

    self->bndNodeSet = NULL;
    self->bndElementSet = NULL;
    self->feElType = elType;
    self->feElFamily = (char*)family;
    self->elementMesh = elementMesh;

    /* checkpoint non-constant meshes */
    if ( self->feElFamily && strcmp( self->feElFamily, "constant" ) ) {
        self->isCheckpointedAndReloaded = True;
    }

    self->inc = IArray_New();
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SphericalFeMesh_Delete( void* feMesh )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;
    /* Delete the parent. */
    _Mesh_Delete( self );
}

void _SphericalFeMesh_Print( void* feMesh, Stream* stream )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    /* Print parent */
    Journal_Printf( stream, "SphericalFeMesh (ptr): (%p)\n", self );
    _Mesh_Print( self, stream );
}

void _SphericalFeMesh_AssignFromXML( void* feMesh, Stg_ComponentFactory* cf, void* data )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    assert( self );

    _Mesh_AssignFromXML( self, cf, data );

    self->useFeAlgorithms = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"UseFeAlgorithms", True );

    _SphericalFeMesh_Init( self, NULL, Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"elementType", "linear"  ),
                           Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"isElementMesh", False )  );
}

void _SphericalFeMesh_Build( void* feMesh, void* data )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    Stream*		stream;
    ElementType*	elType;
    int nNodes, nEls, e_i, n_i, *nodes;
    IArray *inc=IArray_New();

    assert( self );

    stream = Journal_Register( Info_Type, (Name)self->type  );

    _Mesh_Build( self, data );

    if( !self->elementMesh ) {
       nNodes = Mesh_GetDomainSize( self, MT_VERTEX );
       // build the bndNodeSet from the read in innerSet and outerSet
       self->bndNodeSet = IndexSet_New( nNodes );
       IndexSet_Merge_OR( self->bndNodeSet, self->innerSet );
       IndexSet_Merge_OR( self->bndNodeSet, self->outerSet );
       IndexSet_UpdateMembersCount(self->bndNodeSet);

       nEls = Mesh_GetDomainSize( self, MT_VOLUME );
       self->bndElementSet = IndexSet_New( nEls );
       self->innerElSet = IndexSet_New( nEls );
       self->outerElSet = IndexSet_New( nEls );

       // build the bnd element list
       for( e_i=0; e_i < nEls; e_i++ ) {
          FeMesh_GetElementNodes( self, e_i, inc ); 
          nNodes = IArray_GetSize( inc ); 
          nodes = IArray_GetPtr( inc );
          for( n_i=0; n_i<8; n_i++ ) {
             if( IndexSet_IsMember( self->innerSet, nodes[n_i] ) ) {
                IndexSet_Add( self->innerElSet, e_i );
                break; // bnd element identified, next element
             }
             if( IndexSet_IsMember( self->outerSet, nodes[n_i] ) ) {
                IndexSet_Add( self->outerElSet, e_i );
                break; // bnd element identified, next element
             }
          }
       }
      
       // update element sets;
       IndexSet_UpdateMembersCount( self->outerElSet );
       IndexSet_UpdateMembersCount( self->innerElSet );

       // build boundary element set
       IndexSet_Merge_OR( self->bndElementSet, self->innerElSet );
       IndexSet_Merge_OR( self->bndElementSet, self->outerElSet );
       IndexSet_UpdateMembersCount( self->bndElementSet );
    }

   Stg_Class_Delete( inc );

    Stream_Indent( stream );
    Journal_Printf( stream, "Assigning SphericalFeMesh element types...\n" );
    Stream_Indent( stream );

    if( !strcmp( self->feElFamily, "quadratic" ) ) {
        unsigned	nDims;

        nDims = Mesh_GetDimSize( self );
        if( nDims == 3 )
            elType = (ElementType*)Triquadratic_New( "" );
        else if( nDims == 2 )
            elType = (ElementType*)Biquadratic_New( "" );
        else
            abort();
    } else if( !strcmp( self->feElFamily, "linear" ) ) {
        unsigned	nDims;

        nDims = Mesh_GetDimSize( self );
        if( nDims == 3 )
            elType = (ElementType*)TrilinearElementType_New( "" );
        else if( nDims == 2 )
            elType = (ElementType*)BilinearElementType_New( "" );
        else if( nDims == 1 )
            elType = (ElementType*)LinearElementType_New( "" );
        else
            abort();
    } else if( !strcmp( self->feElFamily, "linear-regular" ) ) {
        unsigned	nDims;

        nDims = Mesh_GetDimSize( self );
        if( nDims == 3 )
            elType = (ElementType*)RegularTrilinear_New( "" );
        else if( nDims == 2 )
            elType = (ElementType*)RegularBilinear_New( "" );
        else
            abort();
    } else if( !strcmp( self->feElFamily, "linear-inner" ) ) {
        unsigned	nDims;

        nDims = Mesh_GetDimSize( self );
        if( nDims == 3 )
            elType = (ElementType*)TrilinearInnerElType_New( "" );
        else if( nDims == 2 )
            elType = (ElementType*)BilinearInnerElType_New( "" );
        else
            abort();
    } else if( !strcmp( self->feElFamily, "constant" ) ) {
        elType = (ElementType*)ConstantElementType_New( "" );
    } else if( !strcmp( self->feElFamily, "p1" ) ) {
        elType = (ElementType*)P1_New( "" );
    } else
        abort();
    SphericalFeMesh_SetElementType( self, elType );
    if( self->feElType )
        Stg_Component_Build( self->feElType, data, False );

    /* If we are using an element mesh we need to use SphericalFeMesh algorithms. */
    if( self->elementMesh ) {
        /*
        		Stg_Class_Delete( self->algorithms );
        		self->algorithms = SphericalFeMesh_Algorithms_New( "", NULL );
        		Mesh_Algorithms_SetMesh( self->algorithms, self );
        		Mesh_Algorithms_Update( self->algorithms );
        */
    }

    if( !self->elementMesh && self->useFeAlgorithms ) {
        /* We need to swap to the SphericalFeMesh element type because the
           geometric versions do not produce the same results. */
        Stg_Class_Delete( self->elTypes[0] );
        self->elTypes[0] =
            (Mesh_ElementType*)FeMesh_ElementType_New();
        Mesh_ElementType_SetMesh( self->elTypes[0], self );
        Mesh_ElementType_Update( self->elTypes[0] );
    }

    
   Grid** grid=NULL;
   unsigned bsSize[3] = {1,2,3};
   // make a bullshit elGrid just for the existing _FiniteElementContext_SaveMesh()
   self->elGridId = ExtensionManager_Add( self->info, (Name)"elementGrid", sizeof(Grid*) );
   grid = (Grid** )Mesh_GetExtension(self, Grid**, self->elGridId );
   *grid = Grid_New( );
   Grid_SetNumDims( *grid, 3 );
   Grid_SetSizes( *grid, bsSize );

    Journal_Printf( stream, "... FE element types are '%s',\n", elType->type );
    Journal_Printf( stream, "... done.\n" );
    Stream_UnIndent( stream );
    Stream_UnIndent( stream );
}

void _SphericalFeMesh_Initialise( void* feMesh, void* data )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    assert( self );

    _Mesh_Initialise( self, data );

    if( self->feElType )
        Stg_Component_Initialise( self->feElType, data, False );
}

void _SphericalFeMesh_Execute( void* feMesh, void* data )
{
}

void _SphericalFeMesh_Destroy( void* feMesh, void* data )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    FreeObject( self->outerSet );
    FreeObject( self->innerSet );
    FreeObject( self->innerElSet );
    FreeObject( self->outerElSet );
    FreeObject( self->bndNodeSet );
    FreeObject( self->bndElementSet );

    SphericalFeMesh_Destruct( self );
    Stg_Class_Delete( self->inc );
    _Mesh_Destroy( self, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void SphericalFeMesh_SetElementType( void* feMesh, ElementType* elType )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    assert( self );

    if( self->feElType ) Stg_Class_Delete( self->feElType );
    self->feElType = elType;
}

void SphericalFeMesh_SetElementFamily( void* feMesh, const char* family )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    assert( self );

    self->feElFamily = (char*)family;
}

ElementType* SphericalFeMesh_GetElementType( void* feMesh, unsigned element )
{
    SphericalFeMesh*	self = (SphericalFeMesh*)feMesh;

    assert( self );

    return self->feElType;
}

unsigned SphericalFeMesh_GetNodeLocalSize( void* feMesh )
{
    return Mesh_GetLocalSize( feMesh, MT_VERTEX );
}

unsigned SphericalFeMesh_GetNodeRemoteSize( void* feMesh )
{
    return Mesh_GetRemoteSize( feMesh, MT_VERTEX );
}

unsigned SphericalFeMesh_GetNodeDomainSize( void* feMesh )
{
    return Mesh_GetDomainSize( feMesh, MT_VERTEX );
}

unsigned SphericalFeMesh_GetNodeGlobalSize( void* feMesh )
{
    return Mesh_GetGlobalSize( feMesh, MT_VERTEX );
}

unsigned SphericalFeMesh_GetElementLocalSize( void* feMesh )
{
    return Mesh_GetLocalSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned SphericalFeMesh_GetElementRemoteSize( void* feMesh )
{
    return Mesh_GetRemoteSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned SphericalFeMesh_GetElementDomainSize( void* feMesh )
{
    return Mesh_GetDomainSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned SphericalFeMesh_GetElementGlobalSize( void* feMesh )
{
    return Mesh_GetGlobalSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned SphericalFeMesh_GetElementNodeSize( void* feMesh, unsigned element )
{
    return Mesh_GetIncidenceSize( feMesh, Mesh_GetDimSize( feMesh ), element, MT_VERTEX );
}

unsigned SphericalFeMesh_GetNodeElementSize( void* feMesh, unsigned node )
{
    return Mesh_GetIncidenceSize( feMesh, MT_VERTEX, node, Mesh_GetDimSize( feMesh ) );
}

void SphericalFeMesh_GetElementNodes( void* feMesh, unsigned element, IArray* inc )
{
    Mesh_GetIncidence( feMesh, Mesh_GetDimSize( feMesh ), element, MT_VERTEX, inc );
}

void SphericalFeMesh_GetNodeElements( void* feMesh, unsigned node, IArray* inc )
{
    Mesh_GetIncidence( feMesh, MT_VERTEX, node, Mesh_GetDimSize( feMesh ), inc );
}

unsigned SphericalFeMesh_ElementDomainToGlobal( void* feMesh, unsigned domain )
{
    return Mesh_DomainToGlobal( feMesh, Mesh_GetDimSize( feMesh ), domain );
}

Bool SphericalFeMesh_ElementGlobalToDomain( void* feMesh, unsigned global, unsigned* domain )
{
    return Mesh_GlobalToDomain( feMesh, Mesh_GetDimSize( feMesh ), global, domain );
}

unsigned SphericalFeMesh_NodeDomainToGlobal( void* feMesh, unsigned domain )
{
    return Mesh_DomainToGlobal( feMesh, MT_VERTEX, domain );
}

Bool SphericalFeMesh_NodeGlobalToDomain( void* feMesh, unsigned global, unsigned* domain )
{
    return Mesh_GlobalToDomain( feMesh, MT_VERTEX, global, domain );
}

void SphericalFeMesh_CoordGlobalToLocal( void* feMesh, unsigned element, double* global, double* local )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    ElementType*	elType;

    assert( self );
    assert( element < SphericalFeMesh_GetElementDomainSize( self ) );
    assert( global );
    assert( local );

    elType = SphericalFeMesh_GetElementType( self, element );
    ElementType_ConvertGlobalCoordToElLocal( elType, self, element, global, local );
}

void SphericalFeMesh_CoordLocalToGlobal( void* feMesh, unsigned element, double* local, double* global )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    unsigned	nDims;
    ElementType*	elType;
    double*		basis;
    int	nElNodes, *elNodes;
    double		dimBasis;
    double*		vert;
    unsigned	n_i, d_i;

    assert( self );
    assert( element < SphericalFeMesh_GetElementDomainSize( self ) );
    assert( global );
    assert( local );

    nDims = Mesh_GetDimSize( self );
    elType = SphericalFeMesh_GetElementType( self, element );
    SphericalFeMesh_GetElementNodes( self, element, self->inc );
    nElNodes = IArray_GetSize( self->inc );
    elNodes = IArray_GetPtr( self->inc );
    basis = AllocArray( double, nElNodes );
    ElementType_EvaluateShapeFunctionsAt( elType, local, basis );

    memset( global, 0, nDims * sizeof(double) );
    for( n_i = 0; n_i < nElNodes; n_i++ ) {
        dimBasis = basis[n_i];
        vert = Mesh_GetVertex( self, elNodes[n_i] );
        for( d_i = 0; d_i < nDims; d_i++ )
            global[d_i] += dimBasis * vert[d_i];
    }

    FreeArray( basis );
}

void SphericalFeMesh_EvalBasis( void* feMesh, unsigned element, double* localCoord, double* basis )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    ElementType*	elType;

    assert( self );
    assert( localCoord );

    elType = SphericalFeMesh_GetElementType( self, element );
    ElementType_EvaluateShapeFunctionsAt( elType, localCoord, basis );
}

void SphericalFeMesh_EvalLocalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    ElementType*	elType;

    assert( self );
    assert( localCoord );
    assert( derivs );

    elType = SphericalFeMesh_GetElementType( self, element );
    ElementType_EvaluateShapeFunctionLocalDerivsAt( elType, localCoord, derivs );
}

void SphericalFeMesh_EvalGlobalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs, double* jacDet )
{
    SphericalFeMesh*		self = (SphericalFeMesh*)feMesh;
    unsigned	nDims;
    ElementType*	elType;
    double		jd;

    assert( self );
    assert( localCoord );
    assert( derivs );

    nDims = Mesh_GetDimSize( self );
    elType = SphericalFeMesh_GetElementType( self, element );
    ElementType_ShapeFunctionsGlobalDerivs( elType, self, element, localCoord, nDims,
                                            &jd, derivs );
    if( jacDet )
        *jacDet = jd;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void SphericalFeMesh_Destruct( SphericalFeMesh* self )
{
    Stg_Class_Delete( self->feElType );
    self->feElFamily = NULL;
    /* Disabling the killing of this object from within this
    component as this will be destroyed by the LiveComponentRegister_DestroyAll function 101109 */
    /*KillObject( self->feElType );*/
}


