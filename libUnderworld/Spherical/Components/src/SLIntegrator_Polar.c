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
#include <Spherical/Spherical.h>

#include "types.h"
#include "SLIntegrator_Polar.h"

#include <assert.h>

/** Textual name of this class */
const Type SLIntegrator_Polar_Type = "SLIntegrator_Polar";

SLIntegrator_Polar* _SLIntegrator_Polar_New( SLINTEGRATOR_POLAR_DEFARGS ) {
    SLIntegrator_Polar*		self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SLIntegrator_Polar) );
    /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
    /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
       and so should be set to ZERO in any children of this class. */
    nameAllocationType = NON_GLOBAL;

    self = (SLIntegrator_Polar*) _Stg_Component_New( STG_COMPONENT_PASSARGS );

    /* General info */
    self->variableList = Stg_ObjectList_New();
    self->varStarList  = Stg_ObjectList_New();

    return self;
}

void* _SLIntegrator_Polar_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    SLIntegrator_Polar*	self 	= (SLIntegrator_Polar*)slIntegrator;
    SLIntegrator_Polar*	newSLIntegrator_Polar;
    PtrMap*			map 	= ptrMap;
    Bool			ownMap 	= False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSLIntegrator_Polar = _Stg_Component_Copy( self, dest, deep, nameExt, map );

    if( deep ) {
        if( (newSLIntegrator_Polar->velocityField = PtrMap_Find( map, self->velocityField )) == NULL ) {
            newSLIntegrator_Polar->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
            PtrMap_Append( map, self->velocityField, newSLIntegrator_Polar->velocityField );
        }
    }
    else {
        newSLIntegrator_Polar->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSLIntegrator_Polar;
}


void _SLIntegrator_Polar_Delete( void* slIntegrator ) {
    SLIntegrator_Polar*		self = (SLIntegrator_Polar*)slIntegrator;
    Stg_Class_Delete( self->variableList );
    Stg_Class_Delete( self->varStarList );
}

void _SLIntegrator_Polar_Print( void* slIntegrator, Stream* stream ) {
    SLIntegrator_Polar*		self = (SLIntegrator_Polar*)slIntegrator;

    _Stg_Component_Print( self, stream );

    Journal_PrintPointer( stream, self->velocityField );
}

void* _SLIntegrator_Polar_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SLIntegrator_Polar);
    Type                                                      type = SLIntegrator_Polar_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLIntegrator_Polar_Delete;
    Stg_Class_PrintFunction*                                _print = _SLIntegrator_Polar_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLIntegrator_Polar_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SLIntegrator_Polar_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _SLIntegrator_Polar_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SLIntegrator_Polar_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SLIntegrator_Polar_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLIntegrator_Polar_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLIntegrator_Polar_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL; /* default value NON_GLOBAL */

    return (void*)_SLIntegrator_Polar_New( SLINTEGRATOR_POLAR_PASSARGS );
}

void _SLIntegrator_Polar_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data ) {
    SLIntegrator_Polar*		self 		= (SLIntegrator_Polar*)slIntegrator;
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

    EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), SLIntegrator_Polar_InitSolve, self );

    sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", Energy_SLE, False, NULL );
    if( sle ) {
        /* also set sle to run where required */
        EP_InsertClassHookAfter( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "SLIntegrator_Polar_InitSolve", 
            SystemLinearEquations_GetRunEPFunction(), sle );
        /* remember to disable the standard run at execute */
        SystemLinearEquations_SetRunDuringExecutePhase( sle, False );
    }

    self->isConstructed = True;
}

void _SLIntegrator_Polar_Build( void* slIntegrator, void* data ) {
    SLIntegrator_Polar*	self 		= (SLIntegrator_Polar*)slIntegrator;
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

void _SLIntegrator_Polar_Initialise( void* slIntegrator, void* data ) {
    SLIntegrator_Polar*		self 		= (SLIntegrator_Polar*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i;

    if( self->velocityField ) Stg_Component_Initialise( self->velocityField, data, False );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Initialise( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Initialise( feVarStar , data, False );
    }

    self->abcissa = malloc(4*sizeof(double));
    self->Ni      = malloc(64*sizeof(double));
    self->GNix    = malloc(2*sizeof(double*));
    self->GNix[0] = malloc(64*sizeof(double));
    self->GNix[1] = malloc(64*sizeof(double));

    self->abcissa[0] = -1.0;
    self->abcissa[1] = -0.33333333333333333;
    self->abcissa[2] = -1.0*self->abcissa[1];
    self->abcissa[3] = -1.0*self->abcissa[0];

    /* sanity checks */
    if( Mesh_GetDimSize( self->velocityField->feMesh ) != 2 ) {
        printf( "ERROR: component %s requires mesh of dimension 2\n", self->type );
        abort();
    }
    if( strcmp( self->velocityField->feMesh->algorithms->type, "Mesh_SphericalAlgorithms" ) != 0 ) {
        printf( "ERROR: component %s required mesh algorithms type Mesh_SphericalAlgoritms, type is %s.\n", self->velocityField->feMesh->algorithms->type );
        abort();
    }
}

void _SLIntegrator_Polar_Execute( void* slIntegrator, void* data ) {}

void _SLIntegrator_Polar_Destroy( void* slIntegrator, void* data ) {
    SLIntegrator_Polar*		self 		= (SLIntegrator_Polar*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i;

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
    free(self->GNix);

    Stg_Class_Delete( self->inc );
}

void SLIntegrator_Polar_InitSolve( void* _self, void* _context ) {
    SLIntegrator_Polar*		self			= (SLIntegrator_Polar*) _self;
    unsigned			field_i, node_i;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    double            		phi[3];

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];

        /* generate the _* field */
        SLIntegrator_Polar_Solve( self, feVariable, feVarStar );

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
void SLIntegrator_Polar_IntegrateRK4( void* slIntegrator, FeVariable* velocityField, double dt, double* origin, double* position ) {
    SLIntegrator_Polar*	self 	     = (SLIntegrator_Polar*)slIntegrator;
    unsigned		ndims		     = Mesh_GetDimSize( velocityField->feMesh );
    unsigned		dim_i;
    double		min[3], max[3];
    double		k[4][3];
    double		coordPrime[3];
    unsigned*		periodic	     = ((CartesianGenerator*)velocityField->feMesh->generator)->periodic;

    Mesh_GetGlobalCoordRange( velocityField->feMesh, min, max );

    SLIntegrator_Polar_CubicInterpolator( self, velocityField, origin, k[0] );
    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[0][dim_i];
    }
    //if( !spherical ) SLIntegrator_Polar_PeriodicUpdate( coordPrime, min, max, periodic, ndims );
    SLIntegrator_Polar_CubicInterpolator( self, velocityField, coordPrime, k[1] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[1][dim_i];
    }
    //if( !spherical ) SLIntegrator_Polar_PeriodicUpdate( coordPrime, min, max, periodic, ndims );
    SLIntegrator_Polar_CubicInterpolator( self, velocityField, coordPrime, k[2] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - dt*k[2][dim_i];
    }
    //if( !spherical ) SLIntegrator_Polar_PeriodicUpdate( coordPrime, min, max, periodic, ndims );
    SLIntegrator_Polar_CubicInterpolator( self, velocityField, coordPrime, k[3] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        position[dim_i] = origin[dim_i] - INV6*dt*( k[0][dim_i] + 2.0*k[1][dim_i] + 2.0*k[2][dim_i] + k[3][dim_i] );
    }
    //if( !spherical ) SLIntegrator_Polar_PeriodicUpdate( position, min, max, periodic, ndims );
}

Bool HasSide( FeMesh* feMesh, IArray* inc, unsigned elInd, unsigned* elNodes, unsigned nNodes, unsigned* sideNodes ) {
    unsigned	nInc;

    Mesh_GetIncidence( feMesh, MT_VERTEX, elNodes[sideNodes[0]], 2, inc );
    nInc   = IArray_GetSize( inc );
    Mesh_GetIncidence( feMesh, MT_VERTEX, elNodes[sideNodes[1]], 2, inc );
    nInc  += IArray_GetSize( inc );
    if( nInc > 5 )
        return True;

    return False;
}

Bool HasLeft( FeMesh* feMesh, IArray* inc, unsigned elInd, unsigned* elNodes, unsigned nNodes ) {
    Bool 	quad		= (nNodes%3==0) ? True : False;
    unsigned	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (quad) ? 6 : 2;
    return HasSide( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool HasRight( FeMesh* feMesh, IArray* inc, unsigned elInd, unsigned* elNodes, unsigned nNodes ) {
    Bool 	quad		= (nNodes%3==0) ? True : False;
    unsigned	sideNodes[4];

    sideNodes[0] = (quad) ? 2 : 1;
    sideNodes[1] = (quad) ? 8 : 3;
    return HasSide( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool HasBottom( FeMesh* feMesh, IArray* inc, unsigned elInd, unsigned* elNodes, unsigned nNodes ) {
    Bool 	quad		= (nNodes%3==0) ? True : False;
    unsigned	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (quad) ? 2 : 1;
    return HasSide( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool HasTop( FeMesh* feMesh, IArray* inc, unsigned elInd, unsigned* elNodes, unsigned nNodes ) {
    Bool 	quad		= (nNodes%3==0) ? True : False;
    unsigned	sideNodes[4];

    sideNodes[0] = (quad) ? 6 : 2;
    sideNodes[1] = (quad) ? 8 : 3;
    return HasSide( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SLIntegrator_Polar_PeriodicUpdate( double* pos, double* min, double* max, Bool* isPeriodic ) {
    Bool	result		= False;
    unsigned 	dim_i;

    for( dim_i = 0; dim_i < 2; dim_i++ ) {
        if( pos[dim_i] < min[dim_i] ) {
            pos[dim_i] = (isPeriodic[dim_i]) ? max[dim_i] - min[dim_i] + pos[dim_i] : min[dim_i];
            result = True;
        }
        if( pos[dim_i] > max[dim_i] ) {
            pos[dim_i] = (isPeriodic[dim_i]) ? min[dim_i] - max[dim_i] + pos[dim_i] : max[dim_i];
            result = True;
        }
    }

    return result;
}

void SLIntegrator_Polar_GlobalNodeToIJK( FeMesh* feMesh, unsigned gNode, unsigned* IJK ) {
    Grid**      grid            = (Grid**) Mesh_GetExtension( feMesh, Grid*, feMesh->vertGridId );
    unsigned*   sizes           = Grid_GetSizes( *grid );

    IJK[0] = gNode%sizes[0];
    IJK[1] = (gNode/sizes[0])%sizes[1];
}

unsigned SLIntegrator_Polar_IJKToGlobalNode( FeMesh* feMesh, unsigned* IJK ) {
    Grid**      grid            = (Grid**) Mesh_GetExtension( feMesh, Grid*, feMesh->vertGridId );
    unsigned*   sizes           = Grid_GetSizes( *grid );
    unsigned	gNode;

    gNode = sizes[0]*IJK[1] + IJK[0];

    return gNode;
}

void SLIntegrator_Polar_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* position, double* result ) {
    SLIntegrator_Polar*	self 	= (SLIntegrator_Polar*)slIntegrator;
    FeMesh*	feMesh			= feVariable->feMesh;
    unsigned	elInd, nInc, *inc;
    unsigned	IJK[3], IJK_test[3], index = 0;
    unsigned	gNode_I, lNode_I;
    unsigned	nDims			= Mesh_GetDimSize( feMesh );
    unsigned	nodeIndex[64];
    unsigned	numdofs			= feVariable->dofLayout->dofCounts[0];
    double	lCoord[3], phi_i[3];
    unsigned	node_i, dof_i;
    Grid**      grid            	= (Grid**) Mesh_GetExtension( feMesh, Grid*, feMesh->vertGridId );
    unsigned*   sizes           	= Grid_GetSizes( *grid );

    Mesh_SearchElements( feMesh, position, &elInd );
    FeMesh_GetElementNodes( feMesh, elInd, feVariable->inc );
    nInc = IArray_GetSize( feVariable->inc );
    inc  = IArray_GetPtr( feVariable->inc );

    gNode_I = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
    SLIntegrator_Polar_GlobalNodeToIJK( feMesh, gNode_I, IJK );

    if( nInc%3 == 0 ) { /* quadratic mesh */
        unsigned cIndex = ( nDims == 3 ) ? 13 : 4;
        double*  cCoord = Mesh_GetVertex( feMesh, inc[cIndex] );

        if( position[0] > cCoord[0] && !HasRight( feMesh, self->inc, elInd, inc, nInc ) ) IJK[0] -= 2;
        if( position[1] > cCoord[1] && !HasTop( feMesh, self->inc, elInd, inc, nInc ) )   IJK[1] -= 2;
        if( position[0] < cCoord[0] && HasLeft( feMesh, self->inc, elInd, inc, nInc ) )   IJK[0]--;
        if( position[1] < cCoord[1] && HasBottom( feMesh, self->inc, elInd, inc, nInc ) ) IJK[1]--;
    }
    else { /* linear mesh */
        if( !HasRight( feMesh, self->inc, elInd, inc, nInc ) ) IJK[0]--;
        if( !HasTop( feMesh, self->inc, elInd, inc, nInc ) )   IJK[1]--;
        if( HasLeft( feMesh, self->inc, elInd, inc, nInc ) )   IJK[0]--;
        if( HasBottom( feMesh, self->inc, elInd, inc, nInc ) ) IJK[1]--;
    }
    FeMesh_NodeGlobalToDomain( feMesh, SLIntegrator_Polar_IJKToGlobalNode( feMesh, IJK ), &lNode_I );

    /* interpolate using SLIntegrator_Polar_Lagrange polynomial shape functions */
    for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
        for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
            gNode_I = SLIntegrator_Polar_IJKToGlobalNode( feMesh, IJK_test );
            if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_I, &lNode_I ) ) abort();
            else
                nodeIndex[index++] = lNode_I;
        }
    }
    /* map to master element and interpolate */
    SLIntegrator_Polar_GlobalToLocal( self, feMesh, nodeIndex, position, lCoord );
    SLIntegrator_Polar_ShapeFuncs( self, lCoord, self->Ni );
    for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
        result[dof_i] = 0.0;
    }
    for( node_i = 0; node_i < 16; node_i++ ) {
        FeVariable_GetValueAtNode( feVariable, nodeIndex[node_i], phi_i );
        for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
            result[dof_i] += phi_i[dof_i]*self->Ni[node_i];
        }
    }
}

void SLIntegrator_Polar_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* varStarField ) {
    SLIntegrator_Polar*		self 		     = (SLIntegrator_Polar*)slIntegrator;
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
	coord = Mesh_GetVertex(feMesh,node_I);

        SLIntegrator_Polar_IntegrateRK4( self, velocityField, dt, coord, position );
        SLIntegrator_Polar_CubicInterpolator( self, variableField, position, var );
        FeVariable_SetValueAtNode( varStarField, node_I, var );
    }
    FeVariable_SyncShadowValues( varStarField );
}

double SLIntegrator_Polar_Lagrange( double* xi, double x, int j ) {
    int i;
    double l = 1.0;

    for( i = 0; i < 4; i++ ) {
        if( i == j )
            continue;

        l *= (x - xi[i])/(xi[j]-xi[i]);
    }
    return l;
}

double SLIntegrator_Polar_LagrangeDeriv( double* xi, double x, int j ) {
    int i, k;
    double alpha = 1.0, b = 0.0, c = 0.0;

    for( i = 0; i < 4; i++ ) {
        if( i == j )
            continue;

        alpha *= (xi[j]-xi[i]);
        b     += xi[i];

        for( k = i + 1; k < 4; k++ ) {
            if( k == j )
                continue;

            c += xi[i]*xi[k];
        }
    }
    return (3.0*x*x - 2.0*b*x + c)/alpha;
}

void SLIntegrator_Polar_ShapeFuncs( void* slIntegrator, const double xi[], double* const Ni ) {
    SLIntegrator_Polar* self = (SLIntegrator_Polar*) slIntegrator;
    int                     x_j, y_j, pt_j;
    double                  Nx, Ny;

    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        Nx  = SLIntegrator_Polar_Lagrange( self->abcissa, xi[0], x_j );
        Ny  = SLIntegrator_Polar_Lagrange( self->abcissa, xi[1], y_j );

        Ni[pt_j] = Nx*Ny;
    }
}

void SLIntegrator_Polar_ShapeFuncDerivs( void* slIntegrator, const double xi[], double** const GNix ) {
    SLIntegrator_Polar* self = (SLIntegrator_Polar*) slIntegrator;
    double                  Nx, Ny, GNx, GNy;
    int                     x_j, y_j, pt_j;

    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        Nx  = SLIntegrator_Polar_Lagrange( self->abcissa, xi[0], x_j );
        Ny  = SLIntegrator_Polar_Lagrange( self->abcissa, xi[1], y_j );
        GNx = SLIntegrator_Polar_LagrangeDeriv( self->abcissa, xi[0], x_j );
        GNy = SLIntegrator_Polar_LagrangeDeriv( self->abcissa, xi[1], y_j );

        GNix[0][pt_j] = GNx*Ny;
        GNix[1][pt_j] = Nx*GNy;
    }
}

void SLIntegrator_Polar_GlobalToLocal( void* slIntegrator, void* _mesh, unsigned* nodeInds, const double* gCoord, double* lCoord ) {
    SLIntegrator_Polar* 	self 		= (SLIntegrator_Polar*)slIntegrator;
    Mesh*			mesh 		= (Mesh*)_mesh;
    TensorArray         	jacobiMatrix;
    double              	tolerance       = 0.0001;
    double              	maxResidual;
    Iteration_Index     	maxIterations   = 100;
    Iteration_Index     	iteration_I;
    Node_Index          	node_I;
    Node_Index          	nodeCount       = 16;
    XYZ                 	rightHandSide;
    XYZ                 	xiIncrement;
    double*       	    	nodeCoord;
    double*		    	Ni		= self->Ni;
    double**	    		GNix		= self->GNix;

    /* Initial guess for element local coordinate is in the centre of the element - ( 0.0, 0.0, 0.0 ) */
    memset( lCoord, 0, 2*sizeof(double) );

    /* Do Newton-Raphson Iteration */
    for ( iteration_I = 0 ; iteration_I < maxIterations ; iteration_I++ ) {
        /* Initialise Values */
        TensorArray_Zero( jacobiMatrix );
        memset( rightHandSide, 0, sizeof( XYZ ) );

        /* Evaluate shape functions for rhs */
        SLIntegrator_Polar_ShapeFuncs( self, lCoord, Ni );
        SLIntegrator_Polar_ShapeFuncDerivs( self, lCoord, GNix );

        for( node_I = 0 ; node_I < nodeCount ; node_I++ ) {
            nodeCoord = Mesh_GetVertex( mesh, nodeInds[node_I] );

            /* Form jacobi matrix */
            jacobiMatrix[ MAP_TENSOR( 0, 0, 2 ) ] += GNix[0][node_I] * nodeCoord[ I_AXIS ];
            jacobiMatrix[ MAP_TENSOR( 0, 1, 2 ) ] += GNix[1][node_I] * nodeCoord[ I_AXIS ];

            jacobiMatrix[ MAP_TENSOR( 1, 0, 2 ) ] += GNix[0][node_I] * nodeCoord[ J_AXIS ];
            jacobiMatrix[ MAP_TENSOR( 1, 1, 2 ) ] += GNix[1][node_I] * nodeCoord[ J_AXIS ];

            /* Form right hand side */
            rightHandSide[ I_AXIS ] -= Ni[node_I] * nodeCoord[ I_AXIS ];
            rightHandSide[ J_AXIS ] -= Ni[node_I] * nodeCoord[ J_AXIS ];
        }

        /* Finish building right hand side */
        rightHandSide[ I_AXIS ] += gCoord[ I_AXIS ];
        rightHandSide[ J_AXIS ] += gCoord[ J_AXIS ];

        /* Solve for xi increment */
        TensorArray_SolveSystem( jacobiMatrix, xiIncrement, rightHandSide, 2 );

        /* Update xi */
        lCoord[ I_AXIS ] += xiIncrement[ I_AXIS ];
	lCoord[ J_AXIS ] += xiIncrement[ J_AXIS ];

        /* Check for convergence */
        maxResidual = fabs( xiIncrement[ I_AXIS ] );
        if( maxResidual < fabs( xiIncrement[ J_AXIS ] ) )
            maxResidual = fabs( xiIncrement[ J_AXIS ] );

        if( maxResidual < tolerance )
            return;
	}
        /* if we are here, it means the iterative method didn't converge.
           Thus we set the local coord's to be invalid, i.e. greater than 1.0 */
        lCoord[0] = 1.1;
}
