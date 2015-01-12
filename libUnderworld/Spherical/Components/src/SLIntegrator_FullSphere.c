/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	AuScope - http://www.auscope.org
**	Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "Components.h"

#include <assert.h>

#define DIM 3

/** Textual name of this class */
const Type SLIntegrator_FullSphere_Type = "SLIntegrator_FullSphere";

SLIntegrator_FullSphere* _SLIntegrator_FullSphere_New( SLINTEGRATOR_FULLSPHERE_DEFARGS ) {
    SLIntegrator_FullSphere*		self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SLIntegrator_FullSphere) );
    /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
    /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
       and so should be set to ZERO in any children of this class. */
    nameAllocationType = NON_GLOBAL;

    self = (SLIntegrator_FullSphere*) _Stg_Component_New( STG_COMPONENT_PASSARGS );

    /* General info */
    self->variableList = Stg_ObjectList_New();
    self->varStarList  = Stg_ObjectList_New();

    return self;
}

void* _SLIntegrator_FullSphere_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    SLIntegrator_FullSphere*	self 	= (SLIntegrator_FullSphere*)slIntegrator;
    SLIntegrator_FullSphere*	newSLIntegrator_FullSphere;
    PtrMap*			map 	= ptrMap;
    Bool			ownMap 	= False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSLIntegrator_FullSphere = _Stg_Component_Copy( self, dest, deep, nameExt, map );

    if( deep ) {
        if( (newSLIntegrator_FullSphere->velocityField = PtrMap_Find( map, self->velocityField )) == NULL ) {
            newSLIntegrator_FullSphere->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
            PtrMap_Append( map, self->velocityField, newSLIntegrator_FullSphere->velocityField );
        }
    }
    else {
        newSLIntegrator_FullSphere->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSLIntegrator_FullSphere;
}


void _SLIntegrator_FullSphere_Delete( void* slIntegrator ) {
    SLIntegrator_FullSphere*		self = (SLIntegrator_FullSphere*)slIntegrator;
    Stg_Class_Delete( self->variableList );
    Stg_Class_Delete( self->varStarList );
}

void _SLIntegrator_FullSphere_Print( void* slIntegrator, Stream* stream ) {
    SLIntegrator_FullSphere*		self = (SLIntegrator_FullSphere*)slIntegrator;

    _Stg_Component_Print( self, stream );

    Journal_PrintPointer( stream, self->velocityField );
}

void* _SLIntegrator_FullSphere_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SLIntegrator_FullSphere);
    Type                                                      type = SLIntegrator_FullSphere_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLIntegrator_FullSphere_Delete;
    Stg_Class_PrintFunction*                                _print = _SLIntegrator_FullSphere_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLIntegrator_FullSphere_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SLIntegrator_FullSphere_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _SLIntegrator_FullSphere_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SLIntegrator_FullSphere_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SLIntegrator_FullSphere_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLIntegrator_FullSphere_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLIntegrator_FullSphere_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL; /* default value NON_GLOBAL */

    return (void*)_SLIntegrator_FullSphere_New( SLINTEGRATOR_FULLSPHERE_PASSARGS );
}

void _SLIntegrator_FullSphere_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data ) {
    SLIntegrator_FullSphere*	self 		= (SLIntegrator_FullSphere*)slIntegrator;
    Dictionary*			dict;
    Dictionary_Entry_Value*	dev;
    unsigned			field_i;
    Name			fieldName;
    FeVariable*	   		feVariable;
    unsigned			nEntries;
    Energy_SLE*			sle		= NULL;

    Stg_Component_AssignFromXML( self, cf, data, False );

    self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", FiniteElementContext, False, data );
    if( !self->context  )
        self->context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  );

    self->velocityField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"VelocityField", FeVariable, False, NULL );

    dict     = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, (Dictionary_Entry_Key)self->name )  );
    dev      = Dictionary_Get( dict, (Dictionary_Entry_Key)"fields" );
    nEntries = Dictionary_Entry_Value_GetCount( dev );

    self->pureAdvection = Memory_Alloc_Array_Unnamed( Bool, nEntries/3 );


    for( field_i = 0; field_i < nEntries; field_i += 3  ) {
        fieldName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( dev, field_i ) );
        feVariable = Stg_ComponentFactory_ConstructByName( cf, (Name)fieldName, FeVariable, True, data );
        Stg_ObjectList_Append( self->variableList, feVariable );

        /* the corresponding _* field */
        fieldName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( dev, field_i + 1 ) );
        feVariable = Stg_ComponentFactory_ConstructByName( cf, (Name)fieldName, FeVariable, True, data );
        Stg_ObjectList_Append( self->varStarList, feVariable );

        /* the corresponding boolean for if the field is pure advection (ie not part of an SLE and so updated here) */
        self->pureAdvection[field_i/3] = Dictionary_Entry_Value_AsBool( Dictionary_Entry_Value_GetElement( dev, field_i + 2 ) );
    }

    EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), SLIntegrator_FullSphere_InitSolve, self );

    sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", Energy_SLE, False, NULL );
    if( sle ) {
        /* also set sle to run where required */
        EP_InsertClassHookAfter( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "SLIntegrator_FullSphere_InitSolve", 
            SystemLinearEquations_GetRunEPFunction(), sle );
        /* remember to disable the standard run at execute */
        SystemLinearEquations_SetRunDuringExecutePhase( sle, False );
    }

    self->isConstructed = True;
}

void _SLIntegrator_FullSphere_Build( void* slIntegrator, void* data ) {
    SLIntegrator_FullSphere*	self 		= (SLIntegrator_FullSphere*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i;

    if(self->velocityField) Stg_Component_Build(self->velocityField, data, False);

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if(feVariable) Stg_Component_Build(feVariable, data, False);
        if(feVarStar ) Stg_Component_Build(feVarStar , data, False);
    }

    self->inc = IArray_New();
}

void _SLIntegrator_FullSphere_Initialise( void* slIntegrator, void* data ) {
    SLIntegrator_FullSphere*	self 		= (SLIntegrator_FullSphere*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    unsigned			field_i, el_i;

    if( self->velocityField ) Stg_Component_Initialise( self->velocityField, data, False );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Initialise( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Initialise( feVarStar , data, False );
    }

    self->abcissa = malloc(4*sizeof(double));
    self->Ni      = malloc(64*sizeof(double));
    self->GNix    = malloc(DIM*sizeof(double*));
    self->GNix[0] = malloc(64*sizeof(double));
    self->GNix[1] = malloc(64*sizeof(double));
    self->GNix[2] = malloc(64*sizeof(double));

    self->elPatch = malloc(Mesh_GetLocalSize( feMesh, DIM )*sizeof(unsigned*));
    for( el_i = 0; el_i < Mesh_GetLocalSize( feMesh, DIM ); el_i++ ) {
        self->elPatch[el_i] = malloc( 64*sizeof(unsigned) );
    }

    self->abcissa[0] = -1.0;
    self->abcissa[1] = -0.33333333333333333;
    self->abcissa[2] = -1.0*self->abcissa[1];
    self->abcissa[3] = -1.0*self->abcissa[0];

    /* sanity checks */
    if( Mesh_GetDimSize( feMesh ) != DIM ) {
        printf( "ERROR: component %s requires mesh of dimension 3\n", self->type );
        abort();
    }
    if( strcmp( feMesh->algorithms->type, "Mesh_ProjectionAlgorithms" ) != 0 ) {
        printf( "ERROR: component %s required mesh algorithms type Mesh_SphericalAlgoritms, type is %s.\n", self->type, feMesh->algorithms->type );
        abort();
    }

    SLIntegrator_FullSphere_InitPatches( self );
}

void _SLIntegrator_FullSphere_Execute( void* slIntegrator, void* data ) {}

void _SLIntegrator_FullSphere_Destroy( void* slIntegrator, void* data ) {
    SLIntegrator_FullSphere*	self 		= (SLIntegrator_FullSphere*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i, el_i;

    if( self->velocityField ) Stg_Component_Destroy( self->velocityField, data, False );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Destroy( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Destroy( feVarStar , data, False );
    }
    if( self->pureAdvection ) Memory_Free( self->pureAdvection );

    free(self->abcissa);
    free(self->Ni);
    free(self->GNix[0]);
    free(self->GNix[1]);
    free(self->GNix[2]);
    free(self->GNix);

    for( el_i = 0; el_i < Mesh_GetLocalSize( self->velocityField->feMesh, DIM ); el_i++ ) {
        free( self->elPatch[el_i] );
    }
    free( self->elPatch );

    Stg_Class_Delete( self->inc );
}

void SLIntegrator_FullSphere_InitSolve( void* _self, void* _context ) {
    SLIntegrator_FullSphere*	self			= (SLIntegrator_FullSphere*) _self;
    unsigned			field_i, node_i;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    double            		phi[3];

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];

        /* generate the _* field */
        SLIntegrator_FullSphere_Solve( self, feVariable, feVarStar );

        /* if the field is not part of an SLE, (ie: is purely advected), then update here */
        if( self->pureAdvection[field_i] ) {
            for( node_i = 0; node_i < Mesh_GetLocalSize( feVariable->feMesh, MT_VERTEX ); node_i++ ) {
                FeVariable_GetValueAtNode( feVarStar,  node_i, phi );
                FeVariable_SetValueAtNode( feVariable, node_i, phi );
            }
            FeVariable_SyncShadowValues( feVariable );
        }
    }
}

#define INV6 0.16666666666666667
void SLIntegrator_FullSphere_IntegrateRK4( void* slIntegrator, FeVariable* velocityField, double dt, double* origin, double* position ) {
    SLIntegrator_FullSphere*	self 	     	= (SLIntegrator_FullSphere*)slIntegrator;
    unsigned			dim_i;
    double			k[4][3];
    double			coordPrime[3];

    SLIntegrator_FullSphere_CubicInterpolator( self, velocityField, origin, k[0] );
    for( dim_i = 0; dim_i < DIM; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[0][dim_i];
    }
    SLIntegrator_FullSphere_CubicInterpolator( self, velocityField, coordPrime, k[1] );

    for( dim_i = 0; dim_i < DIM; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[1][dim_i];
    }
    SLIntegrator_FullSphere_CubicInterpolator( self, velocityField, coordPrime, k[2] );

    for( dim_i = 0; dim_i < DIM; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - dt*k[2][dim_i];
    }
    SLIntegrator_FullSphere_CubicInterpolator( self, velocityField, coordPrime, k[3] );

    for( dim_i = 0; dim_i < DIM; dim_i++ ) {
        position[dim_i] = origin[dim_i] - INV6*dt*( k[0][dim_i] + 2.0*k[1][dim_i] + 2.0*k[2][dim_i] + k[3][dim_i] );
    }
}

void SLIntegrator_FullSphere_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* position, double* result ) {
    SLIntegrator_FullSphere*	self 			= (SLIntegrator_FullSphere*)slIntegrator;
    FeMesh*			feMesh			= feVariable->feMesh;
    unsigned			elInd, *nodes;
    unsigned			numdofs			= feVariable->dofLayout->dofCounts[0];
    double			lCoord[3], phi_i[3];
    unsigned			node_i, dof_i;

    Mesh_SearchElements( feMesh, position, &elInd );
    nodes = self->elPatch[elInd];

    SLIntegrator_Spherical_GlobalToLocal( self, feMesh, nodes, position, lCoord );
    SLIntegrator_Spherical_ShapeFuncs3D( self, lCoord, self->Ni );
    for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
        result[dof_i] = 0.0;
    }
    for( node_i = 0; node_i < 64; node_i++ ) {
        FeVariable_GetValueAtNode( feVariable, nodes[node_i], phi_i );
        for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
            result[dof_i] += phi_i[dof_i]*self->Ni[node_i];
        }
    }
}

void SLIntegrator_FullSphere_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* varStarField ) {
    SLIntegrator_FullSphere*	self 		     = (SLIntegrator_FullSphere*)slIntegrator;
    FiniteElementContext*	context		     = self->context;
    unsigned			node_I;
    FeMesh*			feMesh		     = variableField->feMesh;
    unsigned			meshSize	     = Mesh_GetLocalSize( feMesh, MT_VERTEX );
    FeVariable*			velocityField	     = self->velocityField;
    double			dt		     = AbstractContext_Dt( context );
    double			position[3];
    double			var[3];
    double*			coord;

    FeVariable_SyncShadowValues( velocityField );
    FeVariable_SyncShadowValues( variableField );

    /* assume that the variable mesh is the same as the velocity mesh */
    for( node_I = 0; node_I < meshSize; node_I++ ) {
        if( FeVariable_IsBC( variableField, node_I, 0 ) )
            continue;

	coord = Mesh_GetVertex( feMesh, node_I );

        SLIntegrator_FullSphere_IntegrateRK4( self, velocityField, dt, coord, position );
        SLIntegrator_FullSphere_CubicInterpolator( self, variableField, position, var );
        FeVariable_SetValueAtNode( varStarField, node_I, var );
    }
    FeVariable_SyncShadowValues( varStarField );
}

#define SL_LEFT 0
#define SL_RIGHT 1
#define SL_BOTTOM 2
#define SL_TOP 3
#define SL_FRONT 4
#define SL_BACK 5

#define LEFT_BOTTOM_FRONT 0
#define RIGHT_BOTTOM_FRONT 1
#define LEFT_TOP_FRONT 2
#define RIGHT_TOP_FRONT 3
#define LEFT_BOTTOM_BACK 4
#define RIGHT_BOTTOM_BACK 5
#define LEFT_TOP_BACK 6
#define RIGHT_TOP_BACK 7

Bool _HasVertexOnly( unsigned vert, unsigned* elNodes, unsigned* otherElNodes ) {
    unsigned vert_i, vert_j;
    Bool     hasVert	= False;

    for( vert_i = 0; vert_i < 8; vert_i++ ) {
        if( otherElNodes[vert_i] == vert )
            hasVert = True;
    }
    if( !hasVert ) return False;

    for( vert_i = 0; vert_i < 8; vert_i++ ) {
        for( vert_j = 0; vert_j < 8; vert_j++ ) {
            if( elNodes[vert_i] == otherElNodes[vert_j] && elNodes[vert_i] != vert )
                return False;
        }
    }
    return True;
}

Bool _FindSide( unsigned* elemNodes, unsigned* sideNodes, unsigned* side ) {
    unsigned 	node_i;
    Bool	found[4]	= { False, False, False, False };

    for( node_i = 0; node_i < 8; node_i++ ) {
        if( elemNodes[node_i] == sideNodes[side[0]] ) found[0] = True;
        if( elemNodes[node_i] == sideNodes[side[1]] ) found[1] = True;
        if( elemNodes[node_i] == sideNodes[side[2]] ) found[2] = True;
        if( elemNodes[node_i] == sideNodes[side[3]] ) found[3] = True;
    }
    if( found[0] && found[1] && found[2] && found[3] )
        return True;

    return False;
}

int FindSide( unsigned* elemNodes, unsigned* sideNodes ) {
    unsigned		left[4] 	= { 0, 2, 4, 6 };
    unsigned		right[4] 	= { 1, 3, 5, 7 };
    unsigned		bottom[4] 	= { 0, 1, 4, 5 };
    unsigned		top[4] 		= { 2, 3, 6, 7 };
    unsigned		front[4] 	= { 0, 1, 2, 3 };
    unsigned		back[4] 	= { 4, 5, 6, 7 };

    if( _FindSide( elemNodes, sideNodes, left   ) ) return SL_LEFT;
    if( _FindSide( elemNodes, sideNodes, right  ) ) return SL_RIGHT;
    if( _FindSide( elemNodes, sideNodes, bottom ) ) return SL_BOTTOM;
    if( _FindSide( elemNodes, sideNodes, top    ) ) return SL_TOP;
    if( _FindSide( elemNodes, sideNodes, front  ) ) return SL_FRONT;
    if( _FindSide( elemNodes, sideNodes, back   ) ) return SL_BACK;

    return -1;
}

int NeighbourElementVert( FeMesh* feMesh, unsigned el, unsigned vert, unsigned* elNodes ) {
    IArray*	elem 	= IArray_New();
    IArray*	node 	= IArray_New();
    IArray*	node2 	= IArray_New();
    unsigned	*elemInc, nElem, *nodeInc, *node2Inc;
    unsigned	elem_i, node_i;
    int		target	= -1;

    Mesh_GetIncidence( feMesh, MT_VOLUME, el, MT_VOLUME, elem );
    elemInc = IArray_GetPtr( elem );
    nElem   = IArray_GetSize( elem );

    Mesh_GetIncidence( feMesh, MT_VOLUME, el, MT_VERTEX, node );
    nodeInc = IArray_GetPtr( node );

    for( elem_i = 0; elem_i < nElem; elem_i++ ) {
        Mesh_GetIncidence( feMesh, MT_VOLUME, elemInc[elem_i], MT_VERTEX, node2 );
        node2Inc = IArray_GetPtr( node2 );
        if( _HasVertexOnly( nodeInc[vert], nodeInc, node2Inc ) ) {
            for( node_i = 0; node_i < 8; node_i++ ) {
                elNodes[node_i] = node2Inc[node_i];
            }
            target = elemInc[elem_i];
        }
    }

    Stg_Class_Delete( elem );
    Stg_Class_Delete( node );
    Stg_Class_Delete( node2 );

    return target;
}

int NeighbourElementEdge( FeMesh* feMesh, unsigned el, unsigned* edge, unsigned* elNodes ) {
    IArray*	elem		= IArray_New();
    IArray*	node		= IArray_New();
    IArray*	node2		= IArray_New();
    unsigned	nElem, *elemInc, *nodeInc, *node2Inc;
    unsigned	elem_i, node_i, node_j;
    unsigned	front_bottom[2]	= { 0, 1 };
    unsigned	front_top[2]	= { 2, 3 };
    unsigned	back_bottom[2]	= { 4, 5 };
    unsigned	back_top[2]	= { 6, 7 };
    unsigned	front_left[2]	= { 0, 2 };
    unsigned	front_right[2]	= { 1, 3 };
    unsigned	back_left[2]	= { 4, 6 };
    unsigned	back_right[2]	= { 5, 7 };
    unsigned	top_left[2]	= { 2, 6 };
    unsigned	top_right[2]	= { 3, 7 };
    unsigned	bottom_left[2]	= { 0, 4 };
    unsigned	bottom_right[2]	= { 1, 5 };
    Bool	found[2], others;
    int 	target		= -1;

    Mesh_GetIncidence( feMesh, MT_VOLUME, el, MT_VOLUME, elem );
    elemInc = IArray_GetPtr( elem );
    nElem   = IArray_GetSize( elem );

    Mesh_GetIncidence( feMesh, MT_VOLUME, el, MT_VERTEX, node );
    nodeInc = IArray_GetPtr( node );

    for( elem_i = 0; elem_i < nElem; elem_i++ ) {
        Mesh_GetIncidence( feMesh, MT_VOLUME, elemInc[elem_i], MT_VERTEX, node2 );
        node2Inc = IArray_GetPtr( node2 );

        found[0] = found[1] = False;
        for( node_i = 0; node_i < 8; node_i++ ) {
            if( node2Inc[node_i] == nodeInc[edge[0]] )
                found[0] = True;
            if( node2Inc[node_i] == nodeInc[edge[1]] )
                found[1] = True;
        }

        others = False;
        if( found[0] && found[1] ) {
            for( node_i = 0; node_i < 8; node_i++ ) {
                if( node_i == edge[0] )
                    continue;
                    
                if( node_i == edge[1] )
                    continue;
                    
                for( node_j = 0; node_j < 8; node_j++ ) {
                    if( node2Inc[node_j] == nodeInc[node_i] )
                        others = True;
                }
            }

            if( !others ) {
                for( node_i = 0; node_i < 8; node_i++ ) {
                    elNodes[node_i] = node2Inc[node_i];
                }
                target = elemInc[elem_i];
            }
        }
    }

    Stg_Class_Delete( elem );
    Stg_Class_Delete( node );
    Stg_Class_Delete( node2 );

    return target;
}

/* get the neighbouring element and its nodes */
int NeighbourElementFace( FeMesh* feMesh, unsigned el, unsigned* side, unsigned* elNodes ) {
    IArray*	node	= IArray_New();
    IArray*	elem	= IArray_New();
    IArray*	node2	= IArray_New();
    unsigned	*nodeInc, *elemInc, nElems, *node2Inc;
    unsigned	node_i, elem_i, node2_i;
    Bool	hasNodes[4];
    int		target	= -1;

    /* get the incident nodes */
    Mesh_GetIncidence( feMesh, MT_VOLUME, el, MT_VERTEX, node );
    nodeInc = IArray_GetPtr( node );

    for( node_i = 0; node_i < 8; node_i++ ) {
        /* for each node get the incident elements */
        Mesh_GetIncidence( feMesh, MT_VERTEX, nodeInc[node_i], MT_VOLUME, elem );
        elemInc = IArray_GetPtr( elem );
        nElems  = IArray_GetSize( elem );
        for( elem_i = 0; elem_i < nElems; elem_i++ ) {
            /* skip over the existing element */
            if( elemInc[elem_i] == el )
                continue;

            /* for each of these elements get the nodes */
            Mesh_GetIncidence( feMesh, MT_VOLUME, elemInc[elem_i], MT_VERTEX, node2 );
            node2Inc = IArray_GetPtr( node2 );

            hasNodes[0] = hasNodes[1] = hasNodes[2] = hasNodes[3] = False;

            /* set true if the nodes form the face */
            for( node2_i = 0; node2_i < 8; node2_i++ ) {
                if( nodeInc[side[0]] == node2Inc[node2_i] ) hasNodes[0] = True;
                if( nodeInc[side[1]] == node2Inc[node2_i] ) hasNodes[1] = True;
                if( nodeInc[side[2]] == node2Inc[node2_i] ) hasNodes[2] = True;
                if( nodeInc[side[3]] == node2Inc[node2_i] ) hasNodes[3] = True;
            }

            if( hasNodes[0] && hasNodes[1] && hasNodes[2] && hasNodes[3] ) {
                for( node2_i = 0; node2_i < 8; node2_i++ ) {
                    elNodes[node2_i] = node2Inc[node2_i];
                }
                target = elemInc[elem_i];
            }
        }
    }

    Stg_Class_Delete( node );
    Stg_Class_Delete( elem );
    Stg_Class_Delete( node2 );

    /* adjacent element not found, face is a boundary */
    return target;
}

Bool IsBoundaryFace( Bool* isBoundaryNode, unsigned* nodes, unsigned* side ) {
    unsigned 	node_i;

    for( node_i = 0; node_i < 4; node_i++ ) {
        if( !isBoundaryNode[nodes[side[node_i]]] )
            return False;
    }
    return True;
}

Bool IsTriplePtEdge( Bool* isTriplePtNode, unsigned* nodes, unsigned* edge ) {
    if( isTriplePtNode[nodes[edge[0]]] && isTriplePtNode[nodes[edge[1]]] )
        return True;

    return False;
}

unsigned Ind( unsigned z, unsigned y, unsigned x ) {
    return 9*z + 3*y + x;
}

/* construct 4x4x4 node patch for the central element elOrig */
void ConstructPatch( FeMesh* feMesh, int elOrig, int* els, unsigned* patch ) {
    IArray*		iArray		= IArray_New();
    unsigned		*nodes, dummy[8];
    unsigned		xi, yi, zi, offset;
    unsigned		left[4] 	= { 0, 2, 4, 6 };
    unsigned		right[4] 	= { 1, 3, 5, 7 };
    unsigned		bottom[4] 	= { 0, 1, 4, 5 };
    unsigned		top[4] 		= { 2, 3, 6, 7 };
    unsigned		front[4] 	= { 0, 1, 2, 3 };
    unsigned		back[4] 	= { 4, 5, 6, 7 };
    unsigned		front_bottom[2]	= { 0, 1 };
    unsigned		front_top[2]	= { 2, 3 };
    unsigned		back_bottom[2]	= { 4, 5 };
    unsigned		back_top[2]	= { 6, 7 };
    unsigned		front_left[2]	= { 0, 2 };
    unsigned		front_right[2]	= { 1, 3 };
    unsigned		back_left[2]	= { 4, 6 };
    unsigned		back_right[2]	= { 5, 7 };
    unsigned		top_left[2]	= { 2, 6 };
    unsigned		top_right[2]	= { 3, 7 };
    unsigned		bottom_left[2]	= { 0, 4 };
    unsigned		bottom_right[2]	= { 1, 5 };

    /* front */
    els[Ind(0,0,0)] = NeighbourElementVert( feMesh, elOrig, LEFT_BOTTOM_FRONT, dummy );
    els[Ind(0,0,1)] = NeighbourElementEdge( feMesh, elOrig, front_bottom, dummy );
    els[Ind(0,0,2)] = NeighbourElementVert( feMesh, elOrig, RIGHT_BOTTOM_FRONT, dummy );

    els[Ind(0,1,0)] = NeighbourElementEdge( feMesh, elOrig, front_left, dummy );
    els[Ind(0,1,1)] = NeighbourElementFace( feMesh, elOrig, front, dummy );
    els[Ind(0,1,2)] = NeighbourElementEdge( feMesh, elOrig, front_right, dummy );

    els[Ind(0,2,0)] = NeighbourElementVert( feMesh, elOrig, LEFT_TOP_FRONT, dummy );
    els[Ind(0,2,1)] = NeighbourElementEdge( feMesh, elOrig, front_top, dummy );
    els[Ind(0,2,2)] = NeighbourElementVert( feMesh, elOrig, RIGHT_TOP_FRONT, dummy );

    /* center */
    els[Ind(1,0,0)] = NeighbourElementEdge( feMesh, elOrig, bottom_left, dummy );
    els[Ind(1,0,1)] = NeighbourElementFace( feMesh, elOrig, bottom, dummy );
    els[Ind(1,0,2)] = NeighbourElementEdge( feMesh, elOrig, bottom_right, dummy );

    els[Ind(1,1,0)] = NeighbourElementFace( feMesh, elOrig, left, dummy );
    els[Ind(1,1,1)] = elOrig;
    els[Ind(1,1,2)] = NeighbourElementFace( feMesh, elOrig, right, dummy );

    els[Ind(1,2,0)] = NeighbourElementEdge( feMesh, elOrig, top_left, dummy );
    els[Ind(1,2,1)] = NeighbourElementFace( feMesh, elOrig, top, dummy );
    els[Ind(1,2,2)] = NeighbourElementEdge( feMesh, elOrig, top_right, dummy );

    /* back */
    els[Ind(2,0,0)] = NeighbourElementVert( feMesh, elOrig, LEFT_BOTTOM_BACK, dummy );
    els[Ind(2,0,1)] = NeighbourElementEdge( feMesh, elOrig, back_bottom, dummy );
    els[Ind(2,0,2)] = NeighbourElementVert( feMesh, elOrig, RIGHT_BOTTOM_BACK, dummy );

    els[Ind(2,1,0)] = NeighbourElementEdge( feMesh, elOrig, back_left, dummy );
    els[Ind(2,1,1)] = NeighbourElementFace( feMesh, elOrig, back, dummy );
    els[Ind(2,1,2)] = NeighbourElementEdge( feMesh, elOrig, back_right, dummy );

    els[Ind(2,2,0)] = NeighbourElementVert( feMesh, elOrig, LEFT_TOP_BACK, dummy );
    els[Ind(2,2,1)] = NeighbourElementEdge( feMesh, elOrig, back_top, dummy );
    els[Ind(2,2,2)] = NeighbourElementVert( feMesh, elOrig, RIGHT_TOP_BACK, dummy );

    for( zi = 0; zi < 3; zi++ ) {
        for( yi = 0; yi < 3; yi++ ) {
            for( xi = 0; xi < 3; xi++ ) {
                Mesh_GetIncidence( feMesh, 3, els[9*zi+3*yi+xi], 0, iArray );
                nodes  = IArray_GetPtr( iArray );
                offset = xi + 4*yi + 16*zi;
                patch[ 0+offset] = nodes[0];
                patch[ 1+offset] = nodes[1];
                patch[ 4+offset] = nodes[2];
                patch[ 5+offset] = nodes[3];
                patch[16+offset] = nodes[4];
                patch[17+offset] = nodes[5];
                patch[20+offset] = nodes[6];
                patch[21+offset] = nodes[7];
            }
        }
    }
    Stg_Class_Delete( iArray );
}

/* construct a table which contains the 4x4x4 node patches associated with each element for interpolation 
   xyz direction convection is -ve: left,  bottom, front 
                               +ve: right, top,    back  
           6------7 
          /|     /|   y
         / |    / |   ^
        2------3  |   |  z
        |  4---|--5   | /
        | /    | /    |/
        |/     |/     o---->x
        0------1
*/
void SLIntegrator_FullSphere_InitPatches( void* slIntegrator ) {
    SLIntegrator_FullSphere*	self 		= (SLIntegrator_FullSphere*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    unsigned			node_i, elem_i;//, gNode_i;
    IArray*			iArray		= IArray_New();
    unsigned			*inc, nInc, inc2[8];
    Bool*			isBoundaryNode	= malloc( Mesh_GetLocalSize( feMesh, MT_VERTEX )*sizeof(unsigned) );
    Bool*			isTriplePtNode	= malloc( Mesh_GetLocalSize( feMesh, MT_VERTEX )*sizeof(unsigned) );
    int				elOrig, els[27];
    /* face indices */
    unsigned			left[4] 	= { 0, 2, 4, 6 };
    unsigned			right[4] 	= { 1, 3, 5, 7 };
    unsigned			bottom[4] 	= { 0, 1, 4, 5 };
    unsigned			top[4] 		= { 2, 3, 6, 7 };
    unsigned			front[4] 	= { 0, 1, 2, 3 };
    unsigned			back[4] 	= { 4, 5, 6, 7 };
    /* edge indices */
    unsigned			front_bottom[2]	= { 0, 1 };
    unsigned			front_top[2]	= { 2, 3 };
    unsigned			back_bottom[2]	= { 4, 5 };
    unsigned			back_top[2]	= { 6, 7 };
    unsigned			front_left[2]	= { 0, 2 };
    unsigned			front_right[2]	= { 1, 3 };
    unsigned			back_left[2]	= { 4, 6 };
    unsigned			back_right[2]	= { 5, 7 };
    unsigned			top_left[2]	= { 2, 6 };
    unsigned			top_right[2]	= { 3, 7 };
    unsigned			bottom_left[2]	= { 0, 4 };
    unsigned			bottom_right[2]	= { 1, 5 };
    Bool			okMove[6];
    unsigned			sideNodes[4];
    int				side;

    /* determine which nodes are at triple points and boundaries,
       - boundary nodes have 4 incident elements (3 at triple points)
       - internal nodes have 8 incident elements (6 at triple points) */
    for( node_i = 0; node_i < Mesh_GetLocalSize( feMesh, MT_VERTEX ); node_i++ ) {
	//gNode_i = FeMesh_NodeDomainToGlobal( feMesh, node_i );
        Mesh_GetIncidence( feMesh, MT_VERTEX, node_i, MT_VOLUME, iArray ); //TODO: check that this takes the global index and not the local...
        nInc = IArray_GetSize( iArray );
        isBoundaryNode[node_i] = ( nInc < 6 ) ? True : False;
        isTriplePtNode[node_i] = ( nInc % 3 == 0 ) ? True : False;
    }

    for( elem_i = 0; elem_i < Mesh_GetLocalSize( feMesh, MT_VOLUME ); elem_i++ ) {
        Mesh_GetIncidence( feMesh, MT_VOLUME, elem_i, MT_VERTEX, iArray );
        inc = IArray_GetPtr( iArray );
        for( node_i = 0; node_i < 8; node_i++ ) inc2[node_i] = inc[node_i];

        elOrig = elem_i;
        okMove[0] = okMove[1] = okMove[2] = okMove[3] = okMove[4] = okMove[5] = True;

	/* construct trial patch centered on current element - note that boundary faces and 
           triple point edges are orthogonal, so testing for these can be performed independently */
        if( IsBoundaryFace( isBoundaryNode, inc2, right  ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[left[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, left,   inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }
        if( IsBoundaryFace( isBoundaryNode, inc2, left   ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[right[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, right,  inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }
        if( IsBoundaryFace( isBoundaryNode, inc2, top    ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[bottom[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, bottom, inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }
        if( IsBoundaryFace( isBoundaryNode, inc2, bottom ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[top[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, top,    inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }
        if( IsBoundaryFace( isBoundaryNode, inc2, back   ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[front[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, front,  inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }
        if( IsBoundaryFace( isBoundaryNode, inc2, front  ) ) {
            for( node_i = 0; node_i < 4; node_i++ ) sideNodes[node_i] = inc2[back[node_i]];
            elOrig       = NeighbourElementFace( feMesh, elOrig, back,   inc2 );
            side         = FindSide( inc2, sideNodes );
            okMove[side] = False;
        }

        /* move bottom/left/front if triple point edge to the top/right/back */
        if( IsTriplePtEdge( isTriplePtNode, inc2, front_bottom ) ) {
            if( okMove[SL_BACK]   ) elOrig = NeighbourElementFace( feMesh, elOrig, back,   inc2 );
            if( okMove[SL_TOP]    ) elOrig = NeighbourElementFace( feMesh, elOrig, top,    inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, front_top    ) ) {
            if( okMove[SL_BACK]   ) elOrig = NeighbourElementFace( feMesh, elOrig, back,   inc2 );
            if( okMove[SL_BOTTOM] ) elOrig = NeighbourElementFace( feMesh, elOrig, bottom, inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, back_bottom  ) ) {
            if( okMove[SL_FRONT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, front,  inc2 );
            if( okMove[SL_TOP]    ) elOrig = NeighbourElementFace( feMesh, elOrig, top,    inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, back_top     ) ) {
            if( okMove[SL_FRONT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, front,  inc2 );
            if( okMove[SL_BOTTOM] ) elOrig = NeighbourElementFace( feMesh, elOrig, bottom, inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, front_left   ) ) {
            if( okMove[SL_BACK]   ) elOrig = NeighbourElementFace( feMesh, elOrig, back,   inc2 );
            if( okMove[SL_RIGHT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, right,  inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, front_right  ) ) {
            if( okMove[SL_BACK]   ) elOrig = NeighbourElementFace( feMesh, elOrig, back,   inc2 );
            if( okMove[SL_LEFT]   ) elOrig = NeighbourElementFace( feMesh, elOrig, left,   inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, back_left    ) ) {
            if( okMove[SL_FRONT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, front,  inc2 );
            if( okMove[SL_RIGHT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, right,  inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, back_right   ) ) {
            if( okMove[SL_FRONT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, front,  inc2 );
            if( okMove[SL_LEFT]   ) elOrig = NeighbourElementFace( feMesh, elOrig, left,   inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, top_left     ) ) {
            if( okMove[SL_BOTTOM] ) elOrig = NeighbourElementFace( feMesh, elOrig, bottom, inc2 );
            if( okMove[SL_RIGHT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, right,  inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, top_right    ) ) {
            if( okMove[SL_BOTTOM] ) elOrig = NeighbourElementFace( feMesh, elOrig, bottom, inc2 );
            if( okMove[SL_LEFT]   ) elOrig = NeighbourElementFace( feMesh, elOrig, left,   inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, bottom_left  ) ) {
            if( okMove[SL_TOP]    ) elOrig = NeighbourElementFace( feMesh, elOrig, top,    inc2 );
            if( okMove[SL_RIGHT]  ) elOrig = NeighbourElementFace( feMesh, elOrig, right,  inc2 );
        }
        if( IsTriplePtEdge( isTriplePtNode, inc2, bottom_right ) ) {
            if( okMove[SL_TOP]    ) elOrig = NeighbourElementFace( feMesh, elOrig, top,    inc2 );
            if( okMove[SL_LEFT]   ) elOrig = NeighbourElementFace( feMesh, elOrig, left,   inc2 );
        }

        /* finally, no boundaries or triple point edges on the central element, so construct the final version of the patch */
        ConstructPatch( feMesh, elOrig, els, self->elPatch[elem_i] );

        /* sanity check - make sure there are not duplicate elements (means we're not properly accounting for a triple point) */
        unsigned ii, jj;
        Bool elFound = False;
        for( jj = 0; jj < 27; jj++ ) {
            for( ii = 0; ii < 27; ii++ ) {
                if( ii == jj )
                    continue;

                if( els[ii] == els[jj] ) {
                    printf( "\nERROR: element %u interpolation patch construction failure - duplicate elements in patch\n\n", elem_i );
                    abort();
                }
            }

            if( els[jj] == elem_i )
                elFound = True;
        }
        if( !elFound ) {
            printf( "\nERROR: element %u interpolation patch construction failure - principle element not in patch\n\n", elem_i );
            abort();
        }
    }

    free( isBoundaryNode );
    free( isTriplePtNode );

    Stg_Class_Delete( iArray );
}
