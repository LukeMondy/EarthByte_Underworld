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

#include "types.h"
#include "SLIntegrator_Unstructured.h"

#include <assert.h>

/** Textual name of this class */
const Type SLIntegrator_Unstructured_Type = "SLIntegrator_Unstructured";

SLIntegrator_Unstructured* _SLIntegrator_Unstructured_New( SEMILAGRANGIANINTEGRATOR_DEFARGS ) {
    SLIntegrator_Unstructured*		self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SLIntegrator_Unstructured) );
    /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
    /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
       and so should be set to ZERO in any children of this class. */
    nameAllocationType = NON_GLOBAL;

    self = (SLIntegrator_Unstructured*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

    /* General info */
    self->variableList = Stg_ObjectList_New();
    self->varStarList  = Stg_ObjectList_New();

    return self;
}

void* _SLIntegrator_Unstructured_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    SLIntegrator_Unstructured*	self 	= (SLIntegrator_Unstructured*)slIntegrator;
    SLIntegrator_Unstructured*	newSLIntegrator_Unstructured;
    PtrMap*			map 	= ptrMap;
    Bool			ownMap 	= False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSLIntegrator_Unstructured = _Stg_Component_Copy( self, dest, deep, nameExt, map );

    if( deep ) {
        if( (newSLIntegrator_Unstructured->velocityField = PtrMap_Find( map, self->velocityField )) == NULL ) {
            newSLIntegrator_Unstructured->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
            PtrMap_Append( map, self->velocityField, newSLIntegrator_Unstructured->velocityField );
        }
    }
    else {
        newSLIntegrator_Unstructured->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSLIntegrator_Unstructured;
}


void _SLIntegrator_Unstructured_Delete( void* slIntegrator ) {
    SLIntegrator_Unstructured*		self = (SLIntegrator_Unstructured*)slIntegrator;
    Stg_Class_Delete( self->variableList );
    Stg_Class_Delete( self->varStarList );
}

void _SLIntegrator_Unstructured_Print( void* slIntegrator, Stream* stream ) {
    SLIntegrator_Unstructured*		self = (SLIntegrator_Unstructured*)slIntegrator;

    _Stg_Component_Print( self, stream );

    Journal_PrintPointer( stream, self->velocityField );
}

void* _SLIntegrator_Unstructured_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SLIntegrator_Unstructured);
    Type                                                      type = SLIntegrator_Unstructured_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLIntegrator_Unstructured_Delete;
    Stg_Class_PrintFunction*                                _print = _SLIntegrator_Unstructured_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLIntegrator_Unstructured_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SLIntegrator_Unstructured_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _SLIntegrator_Unstructured_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SLIntegrator_Unstructured_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SLIntegrator_Unstructured_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLIntegrator_Unstructured_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLIntegrator_Unstructured_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL; /* default value NON_GLOBAL */

    return (void*)_SLIntegrator_Unstructured_New( SEMILAGRANGIANINTEGRATOR_PASSARGS );
}

void _SLIntegrator_Unstructured_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data ) {
    SLIntegrator_Unstructured*	self 		= (SLIntegrator_Unstructured*)slIntegrator;
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

    EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), SLIntegrator_Unstructured_InitSolve, self );

    sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", Energy_SLE, False, NULL );
    if( sle ) {
        /* also set sle to run where required */
        EP_InsertClassHookAfter( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "SLIntegrator_Unstructured_InitSolve", 
            SystemLinearEquations_GetRunEPFunction(), sle );
        /* remember to disable the standard run at execute */
        SystemLinearEquations_SetRunDuringExecutePhase( sle, False );
    }

    self->isConstructed = True;
}

void _SLIntegrator_Unstructured_Build( void* slIntegrator, void* data ) {
    SLIntegrator_Unstructured*	self 		= (SLIntegrator_Unstructured*)slIntegrator;
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
}

double P4  ( double x )  { return 0.125*(35.0*x*x*x*x - 30.0*x*x + 3.0); }
double dP4 ( double x )  { return  0.25*(70.0*x*x*x - 30.0*x);           }
double d2P4( double x )  { return  0.25*(210.0*x*x - 30.0);              }

void _SLIntegrator_Unstructured_Initialise( void* slIntegrator, void* data ) {
    SLIntegrator_Unstructured*	self 		= (SLIntegrator_Unstructured*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i, i, j;

    if(self->velocityField) Stg_Component_Initialise(self->velocityField, data, False);

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if(feVariable) Stg_Component_Initialise(feVariable, data, False);
        if(feVarStar ) Stg_Component_Initialise(feVarStar , data, False);
    }

    //self->weights = malloc(4*sizeof(double));
    self->abcissa = malloc(4*sizeof(double));
    self->Ni      = malloc(16*sizeof(double));
    self->GNix    = malloc(2*sizeof(double*));
    self->GNix[0] = malloc(16*sizeof(double));
    self->GNix[1] = malloc(16*sizeof(double));
    //self->dCij    = malloc(16*sizeof(double));

    //self->weights[0] = 0.16666666666666667;
    //self->weights[1] = 0.83333333333333333;
    //self->weights[2] = self->weights[1];
    //self->weights[3] = self->weights[0];

    //self->abcissa[0] = -1.0;
    //self->abcissa[1] = -0.44721359549995793;
    //self->abcissa[2] = -1.0*self->abcissa[1];
    //self->abcissa[3] = -1.0*self->abcissa[0];

    self->abcissa[0] = -1.0;
    self->abcissa[1] = -0.33333333333333333;
    self->abcissa[2] = -1.0*self->abcissa[1];
    self->abcissa[3] = -1.0*self->abcissa[0];

    /*for( i = 0; i <= 3; i++ ) {
        for( j = 0; j <= 3; j++ ) {
            if     ( i == j && i == 0 ) { self->dCij[i*4+j] = -0.25*3*4; }
            else if( i == j && i == 3 ) { self->dCij[i*4+j] = +0.25*3*4; }
            else if( i == j )           { self->dCij[i*4+j] = 0.0; }
            else {
                self->dCij[i*4+j] = P4( self->abcissa[i] )/( P4( self->abcissa[j] )*( self->abcissa[i] - self->abcissa[j] ) );
            }
        }
    }*/
}

void _SLIntegrator_Unstructured_Execute( void* slIntegrator, void* data ) {
}

void _SLIntegrator_Unstructured_Destroy( void* slIntegrator, void* data ) {
    SLIntegrator_Unstructured*	self 		= (SLIntegrator_Unstructured*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i;

    if(self->velocityField) Stg_Component_Destroy(self->velocityField, data, False);

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if(feVariable) Stg_Component_Destroy(feVariable, data, False);
        if(feVarStar ) Stg_Component_Destroy(feVarStar , data, False);
    }
    if(self->pureAdvection) Memory_Free(self->pureAdvection);

    //free(self->weights);
    free(self->abcissa);
    free(self->Ni);
    free(self->GNix[0]);
    free(self->GNix[1]);
    free(self->GNix);
    //free(self->dCij);
}

void SLIntegrator_Unstructured_InitSolve( void* _self, void* _context ) {
    SLIntegrator_Unstructured*	self			= (SLIntegrator_Unstructured*) _self;
    unsigned			field_i, node_i;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    double            		phi[3];

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];

        /* generate the _* field */
        SLIntegrator_Unstructured_Solve( self, feVariable, feVarStar );

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
void SLIntegrator_Unstructured_IntegrateRK4( void* slIntegrator, FeVariable* velocityField, double dt, double* delta, unsigned* nnodes, double* origin, double* position ) {
    SLIntegrator_Unstructured*	self 	     = (SLIntegrator_Unstructured*)slIntegrator;
    unsigned		ndims		     = Mesh_GetDimSize( velocityField->feMesh );
    unsigned		dim_i;
    double		min[3], max[3];
    double		k[4][3];
    double		coordPrime[3];
    unsigned*		periodic	     = ((CartesianGenerator*)velocityField->feMesh->generator)->periodic;

    Mesh_GetGlobalCoordRange( velocityField->feMesh, min, max );

    SLIntegrator_Unstructured_CubicInterpolator( self, velocityField, origin, delta, nnodes, k[0] );
    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5 * dt * k[0][dim_i];
    }
    if( SLIntegrator_Unstructured_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SLIntegrator_Unstructured_CubicInterpolator( self, velocityField, coordPrime, delta, nnodes, k[1] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5 * dt * k[1][dim_i];
    }
    if( SLIntegrator_Unstructured_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SLIntegrator_Unstructured_CubicInterpolator( self, velocityField, coordPrime, delta, nnodes, k[2] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - dt * k[2][dim_i];
    }
    if( SLIntegrator_Unstructured_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SLIntegrator_Unstructured_CubicInterpolator( self, velocityField, coordPrime, delta, nnodes, k[3] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        position[dim_i] = origin[dim_i] -
            INV6 * dt * ( k[0][dim_i] + 2.0 * k[1][dim_i] + 2.0 * k[2][dim_i] + k[3][dim_i] );
    }
    if( SLIntegrator_Unstructured_PeriodicUpdate( position, min, max, periodic, ndims ) );
}

Bool SLIntegrator_Unstructured_PeriodicUpdate( double* pos, double* min, double* max, Bool* isPeriodic, unsigned ndim ) {
    Bool	result		= False;
    unsigned 	dim_i;

    for( dim_i = 0; dim_i < ndim; dim_i++ ) {
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

void SLIntegrator_Unstructured_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* position, double* delta, unsigned* nNodes, double* result ) {
    SLIntegrator_Unstructured*	self 		     = (SLIntegrator_Unstructured*)slIntegrator;
    FeMesh*	feMesh			= feVariable->feMesh;
    Index	elementIndex;
    unsigned	nInc, *inc;
    int		x_0, y_0, z_0;
    int		x_i, y_i, z_i;
    double	localMin[3], localMax[3];
    Index	gNode_I, lNode_I;
    double	px[4], py[4], pz[4];
    unsigned	nDims			= Mesh_GetDimSize( feMesh );
    unsigned	nodeIndex[16];
    unsigned	node_I3D[4][4][4];
    unsigned	numdofs			= feVariable->dofLayout->dofCounts[0];

    Mesh_SearchElements( feMesh, position, &elementIndex );
    FeMesh_GetElementNodes( feMesh, elementIndex, feVariable->inc );
    nInc = IArray_GetSize( feVariable->inc );
    inc  = IArray_GetPtr( feVariable->inc );

    Mesh_GetLocalCoordRange( feMesh, localMin, localMax );
    gNode_I = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );

    x_0 = (int) gNode_I % nNodes[0];
    y_0 = ((int) gNode_I / nNodes[0]) % nNodes[1];
    if( nDims == 3 )
        z_0 = (int) gNode_I / ( nNodes[0] * nNodes[1] );

    /* get the number of nodes across and up that the point lies... */
    if( nInc % 3 == 0 ) { /* quadratic mesh */
        /* bottom left corner of stencil is closer to LHS of element */
        if( position[0] <= Mesh_GetVertex( feMesh, inc[1] )[0] ) x_0--;
        /* bottom left corner of stencil is closer to bottom of element */
        if( position[1] <= Mesh_GetVertex( feMesh, inc[3] )[1] ) y_0--;

        /* LHS node is global domain boundary */
        if( position[0] <= localMin[0] + delta[0] ) x_0++;
        /* RHS node is global domain boundary */
        else if( position[0] >= localMax[0] - delta[0] ) x_0--;

        /* top node is global domain boundary */
        if( position[1] <= localMin[1] + delta[1] ) y_0++;
        /* bottom node is global domain boundary */
        else if( position[1] >= localMax[1] - delta[1] ) y_0--;

        if( nDims == 3 ) {
            if( position[2] <= Mesh_GetVertex( feMesh, inc[9] )[2] ) z_0--;

            if( position[2] <= localMin[2] + delta[2] ) z_0++;
            else if( position[2] >= localMax[2] - delta[2] ) z_0--;
        }
    }
    else if ( nInc % 2 == 0 ) { /* linear mesh */
        if( position[0] > localMin[0] + delta[0] ) x_0--;
        if( position[0] >= localMax[0] - delta[0] ) x_0--;

        if( position[1] > localMin[1] + delta[1] ) y_0--;
        if( position[1] >= localMax[1] - delta[1] ) y_0--;

        if( nDims == 3 ) {
            if( position[2] > localMin[2] + delta[2] ) z_0--;
            if( position[2] >= localMax[2] - delta[2] ) z_0--;
        }
    }
    else abort();

    /* interpolate using Lagrange's formula */
    if( nDims == 2 ) {
        for( y_i = 0; y_i < 4; y_i++ ) {
            for( x_i = 0; x_i < 4; x_i++ ) {
                gNode_I = x_0 + x_i + ( y_0 + y_i ) * nNodes[0];
                if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_I, &lNode_I ) ) abort();
                else
                    nodeIndex[y_i*4+x_i] = lNode_I;
            }
        }
        /* map to master element and interpolate */
        double		lCoord[3], phi_i[3];
        unsigned	node_i, dof_i;
        *result = 0.0;
        SLIntegrator_Unstructured_GlobalToLocal( self, feMesh, nodeIndex, position, lCoord );
        SLIntegrator_Unstructured_ShapeFuncs( self, lCoord, self->Ni );
        for( node_i = 0; node_i < 16; node_i++ ) {
            FeVariable_GetValueAtNode( feVariable, nodeIndex[node_i], phi_i );
            for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
                result[dof_i] += phi_i[dof_i]*self->Ni[node_i];
            }
        }
    }
    else {
        for( z_i = 0; z_i < 4; z_i++ ) {
            for( y_i = 0; y_i < 4; y_i++ ) {
                for( x_i = 0; x_i < 4; x_i++ ) {
                    gNode_I = x_0 + x_i + ( y_0 + y_i ) * nNodes[0] + ( z_0 + z_i ) * nNodes[0] * nNodes[1];
                    if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_I, &lNode_I ) ) abort();
                    else
                        node_I3D[x_i][y_i][z_i] = lNode_I;
                }
            }
        }
    }
}

void SLIntegrator_Unstructured_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* varStarField ) {
    SLIntegrator_Unstructured*	self 		     = (SLIntegrator_Unstructured*)slIntegrator;
    FiniteElementContext*	context		     = self->context;
    unsigned			node_I;
    FeMesh*			feMesh		     = variableField->feMesh;
    unsigned			meshSize	     = Mesh_GetLocalSize( feMesh, MT_VERTEX );
    FeVariable*			velocityField	     = self->velocityField;
    double			dt		     = AbstractContext_Dt( context );
    double			position[3];
    double			var[3];
    unsigned			nNodes[3];
    double			delta[3];
    double*			coord;

    SLIntegrator_Unstructured_GetMaxDelta( variableField, delta, nNodes );

    FeVariable_SyncShadowValues( velocityField );
    FeVariable_SyncShadowValues( variableField );

    /* assume that the variable mesh is the same as the velocity mesh */
    for( node_I = 0; node_I < meshSize; node_I++ ) {
	coord = Mesh_GetVertex(feMesh,node_I);

        SLIntegrator_Unstructured_IntegrateRK4( self, velocityField, dt, delta, nNodes, coord, position );
        SLIntegrator_Unstructured_CubicInterpolator( self, variableField, position, delta, nNodes, var );
        FeVariable_SetValueAtNode( varStarField, node_I, var );
    }
    FeVariable_SyncShadowValues( varStarField );
}

void SLIntegrator_Unstructured_GetMaxDelta( FeVariable* variableField, double* delta, unsigned* nNodes ) {
    unsigned            nInc;
    unsigned*           inc;
    unsigned            dim_i;
    FeMesh*             feMesh          = variableField->feMesh;
    unsigned            nDims           = Mesh_GetDimSize( feMesh );
    Grid**              grid            = (Grid**) Mesh_GetExtension( feMesh, Grid*,  feMesh->elGridId );
    unsigned*           sizes           = Grid_GetSizes( *grid );
    unsigned            el_i;
    double              testDelta[3];

    for( dim_i = 0; dim_i < nDims; dim_i++ ) {
        delta[dim_i] = 0.0;
    }
    for( el_i = 0; el_i < Mesh_GetLocalSize( feMesh, nDims ); el_i++ ) {
        FeMesh_GetElementNodes( variableField->feMesh, el_i, variableField->inc );
        nInc = IArray_GetSize( variableField->inc );
        inc  = IArray_GetPtr( variableField->inc );
        if( nInc%3==0 ) { /* quadratic elements */
            testDelta[0] = Mesh_GetVertex( feMesh, inc[1] )[0] - Mesh_GetVertex( feMesh, inc[0] )[0];
            testDelta[1] = Mesh_GetVertex( feMesh, inc[3] )[1] - Mesh_GetVertex( feMesh, inc[0] )[1];
            if( nDims == 3 ) {
                testDelta[2] = Mesh_GetVertex( feMesh, inc[9] )[2] - Mesh_GetVertex( feMesh, inc[0] )[2];
            }
        }
        else {
            testDelta[0] = Mesh_GetVertex( feMesh, inc[1] )[0] - Mesh_GetVertex( feMesh, inc[0] )[0];
            testDelta[1] = Mesh_GetVertex( feMesh, inc[2] )[1] - Mesh_GetVertex( feMesh, inc[0] )[1];
            if( nDims == 3 ) {
                testDelta[2] = Mesh_GetVertex( feMesh, inc[4] )[2] - Mesh_GetVertex( feMesh, inc[0] )[2];
            }
        }
        for( dim_i = 0; dim_i < nDims; dim_i++ ) {
            if( delta[dim_i] < testDelta[dim_i] ) delta[dim_i] = testDelta[dim_i];
        }
    }

    if( nInc%3 == 0 ) { /* quadratic elements */
        for( dim_i = 0; dim_i < nDims; dim_i++ )
            nNodes[dim_i] = 2*sizes[dim_i] + 1;
    }
    else {
        for( dim_i = 0; dim_i < nDims; dim_i++ )
            nNodes[dim_i] = sizes[dim_i] + 1;
    }
}

#define EPS 1.0e-8
Bool IsRoot( double* abcissa, double x, int* i ) {
    int k;

    for( k = 0; k <= 3; k++ ) {
        if( fabs( x - abcissa[k] ) < EPS ) {
            *i = k;
            return True;
        }
    }
    *i = -1;
    return False;
}

double Lagrange( double* xi, double x, int j ) {
    int i;
    double l = 1.0;

    for( i = 0; i < 4; i++ ) {
        if( i == j )
            continue;

        l *= (x - xi[i])/(xi[j]-xi[i]);
    }
    return l;
}

double LagrangeDeriv( double* xi, double x, int j ) {
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

void SLIntegrator_Unstructured_ShapeFuncs( void* slIntegrator, const double xi[], double* const Ni ) {
    SLIntegrator_Unstructured* self = (SLIntegrator_Unstructured*) slIntegrator;
    int                     x_i, y_i, x_j, y_j, pt_j, N = 3;
    double                  Nx, Ny;

/*
    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        if( IsRoot( self->abcissa, xi[0], &x_i ) ) {
            Nx = ( x_i == x_j ) ? 1.0 : 0.0;
        }
        else {
            Nx = -((1.0 - xi[0]*xi[0])/(N*(N + 1)*P4(self->abcissa[x_j])*(xi[0] - self->abcissa[x_j])))*dP4(xi[0]);
        }

        if( IsRoot( self->abcissa, xi[1], &y_i ) ) {
            Ny = ( y_i == y_j ) ? 1.0 : 0.0;
        }
        else {
            Ny = -((1.0 - xi[1]*xi[1])/(N*(N + 1)*P4(self->abcissa[y_j])*(xi[1] - self->abcissa[y_j])))*dP4(xi[1]);
        }
        Ni[pt_j] = Nx*Ny;
    }
*/
    
    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        Nx  = Lagrange( self->abcissa, xi[0], x_j );
        Ny  = Lagrange( self->abcissa, xi[1], y_j );

        Ni[pt_j] = Nx*Ny;
    }
}

void SLIntegrator_Unstructured_ShapeFuncDerivs( void* slIntegrator, const double xi[], double** const GNix ) {
    SLIntegrator_Unstructured* self = (SLIntegrator_Unstructured*) slIntegrator;
    double                  Nx, Ny, GNx, GNy, a, b, c;
    int                     x_i, y_i, x_j, y_j, pt_j;
    unsigned                N = 3;

/*
    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        if( IsRoot( self->abcissa, xi[0], &x_i ) ) {
            Nx  = ( x_i == x_j ) ? 1.0 : 0.0;
            GNx = self->dCij[x_i*(N+1)+x_j];
        }
        else {
            Nx  = -((1.0 - xi[0]*xi[0])/(N*(N + 1)*P4(self->abcissa[x_j])*(xi[0] - self->abcissa[x_j])))*dP4(xi[0]);
            a   = xi[0]*xi[0] + 1.0 - 2.0*xi[0]*self->abcissa[x_j];
            b   = (1.0 - xi[0]*xi[0])*(xi[0] - self->abcissa[x_j]);
            c   = N*(N + 1)*(xi[0] - self->abcissa[x_j])*(xi[0] - self->abcissa[x_j])*P4(self->abcissa[x_j]);
            GNx = (a*dP4(xi[0]) - b*d2P4(xi[0]))/c;
        }

        if( IsRoot( self->abcissa, xi[1], &y_i ) ) {
            Ny  = ( y_i == y_j ) ? 1.0 : 0.0;
            GNy = self->dCij[y_i*(N+1)+y_j];
        }
        else {
            Ny  = -((1.0 - xi[1]*xi[1])/(N*(N + 1)*P4(self->abcissa[y_j])*(xi[1] - self->abcissa[y_j])))*dP4(xi[1]);
            a   = xi[1]*xi[1] + 1.0 - 2.0*xi[1]*self->abcissa[y_j];
            b   = (1.0 - xi[1]*xi[1])*(xi[1] - self->abcissa[y_j]);
            c   = N*(N + 1)*(xi[1] - self->abcissa[y_j])*(xi[1] - self->abcissa[y_j])*P4(self->abcissa[y_j]);
            GNy = (a*dP4(xi[1]) - b*d2P4(xi[1]))/c;
        }

        GNix[0][pt_j] = GNx*Ny;
        GNix[1][pt_j] = Nx*GNy;
    }
*/
    for( pt_j = 0; pt_j < 16; pt_j++ ) {
        x_j = pt_j%4;
        y_j = pt_j/4;

        Nx  = Lagrange( self->abcissa, xi[0], x_j );
        Ny  = Lagrange( self->abcissa, xi[1], y_j );
        GNx = LagrangeDeriv( self->abcissa, xi[0], x_j );
        GNy = LagrangeDeriv( self->abcissa, xi[1], y_j );

        GNix[0][pt_j] = GNx*Ny;
        GNix[1][pt_j] = Nx*GNy;
    }
}

void SLIntegrator_Unstructured_GlobalToLocal( void* slIntegrator, void* _mesh, unsigned* nodeInds, const double* gCoord, double* lCoord ) {
    SLIntegrator_Unstructured* 	self 		= (SLIntegrator_Unstructured*)slIntegrator;
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
    Dimension_Index     	dim             = 2;
    double*		    	Ni		= self->Ni;
    double**	    		GNix		= self->GNix;

    /* Initial guess for element local coordinate is in the centre of the element - ( 0.0, 0.0, 0.0 ) */
    memset( lCoord, 0, dim*sizeof(double) );

    /* Do Newton-Raphson Iteration */
    for ( iteration_I = 0 ; iteration_I < maxIterations ; iteration_I++ ) {
        /* Initialise Values */
        TensorArray_Zero( jacobiMatrix );
        memset( rightHandSide, 0, sizeof( XYZ ) );

        /* Evaluate shape functions for rhs */
        SLIntegrator_Unstructured_ShapeFuncs( self, lCoord, Ni );
        SLIntegrator_Unstructured_ShapeFuncDerivs( self, lCoord, GNix );

        for( node_I = 0 ; node_I < nodeCount ; node_I++ ) {
            nodeCoord = Mesh_GetVertex( mesh, nodeInds[node_I] );

            /* Form jacobi matrix */
            jacobiMatrix[ MAP_TENSOR( 0, 0, dim ) ] += GNix[0][node_I] * nodeCoord[ I_AXIS ];
            jacobiMatrix[ MAP_TENSOR( 0, 1, dim ) ] += GNix[1][node_I] * nodeCoord[ I_AXIS ];

            jacobiMatrix[ MAP_TENSOR( 1, 0, dim ) ] += GNix[0][node_I] * nodeCoord[ J_AXIS ];
            jacobiMatrix[ MAP_TENSOR( 1, 1, dim ) ] += GNix[1][node_I] * nodeCoord[ J_AXIS ];

            /* Form right hand side */
            rightHandSide[ I_AXIS ] -= Ni[node_I] * nodeCoord[ I_AXIS ];
            rightHandSide[ J_AXIS ] -= Ni[node_I] * nodeCoord[ J_AXIS ];
        }

        /* Finish building right hand side */
        rightHandSide[ I_AXIS ] += gCoord[ I_AXIS ];
        rightHandSide[ J_AXIS ] += gCoord[ J_AXIS ];

        /* Solve for xi increment */
        TensorArray_SolveSystem( jacobiMatrix, xiIncrement, rightHandSide, dim );

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
