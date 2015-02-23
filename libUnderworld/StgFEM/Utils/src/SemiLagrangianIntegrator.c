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
#include "SemiLagrangianIntegrator.h"

#include <assert.h>

/** Textual name of this class */
const Type SemiLagrangianIntegrator_Type = "SemiLagrangianIntegrator";

SemiLagrangianIntegrator* _SemiLagrangianIntegrator_New( SEMILAGRANGIANINTEGRATOR_DEFARGS ) {
    SemiLagrangianIntegrator*		self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SemiLagrangianIntegrator) );
    /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
    /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
       and so should be set to ZERO in any children of this class. */
    nameAllocationType = NON_GLOBAL;

    self = (SemiLagrangianIntegrator*) _Stg_Component_New( STG_COMPONENT_PASSARGS );

    /* General info */
    self->variableList = Stg_ObjectList_New();
    self->varStarList  = Stg_ObjectList_New();

    return self;
}

void* _SemiLagrangianIntegrator_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    SemiLagrangianIntegrator*	self 	= (SemiLagrangianIntegrator*)slIntegrator;
    SemiLagrangianIntegrator*	newSemiLagrangianIntegrator;
    PtrMap*			map 	= ptrMap;
    Bool			ownMap 	= False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSemiLagrangianIntegrator = _Stg_Component_Copy( self, dest, deep, nameExt, map );

    if( deep ) {
        if( (newSemiLagrangianIntegrator->velocityField = PtrMap_Find( map, self->velocityField )) == NULL ) {
            newSemiLagrangianIntegrator->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
            PtrMap_Append( map, self->velocityField, newSemiLagrangianIntegrator->velocityField );
        }
    }
    else {
        newSemiLagrangianIntegrator->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSemiLagrangianIntegrator;
}


void _SemiLagrangianIntegrator_Delete( void* slIntegrator ) {
    SemiLagrangianIntegrator*		self = (SemiLagrangianIntegrator*)slIntegrator;
    Stg_Class_Delete( self->variableList );
    Stg_Class_Delete( self->varStarList );
}

void _SemiLagrangianIntegrator_Print( void* slIntegrator, Stream* stream ) {
    SemiLagrangianIntegrator*		self = (SemiLagrangianIntegrator*)slIntegrator;

    _Stg_Component_Print( self, stream );

    Journal_PrintPointer( stream, self->velocityField );
}

void* _SemiLagrangianIntegrator_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SemiLagrangianIntegrator);
    Type                                                      type = SemiLagrangianIntegrator_Type;
    Stg_Class_DeleteFunction*                              _delete = _SemiLagrangianIntegrator_Delete;
    Stg_Class_PrintFunction*                                _print = _SemiLagrangianIntegrator_Print;
    Stg_Class_CopyFunction*                                  _copy = _SemiLagrangianIntegrator_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SemiLagrangianIntegrator_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _SemiLagrangianIntegrator_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SemiLagrangianIntegrator_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SemiLagrangianIntegrator_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SemiLagrangianIntegrator_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SemiLagrangianIntegrator_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL; /* default value NON_GLOBAL */

    return (void*)_SemiLagrangianIntegrator_New( SEMILAGRANGIANINTEGRATOR_PASSARGS );
}

void _SemiLagrangianIntegrator_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
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

    EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), SemiLagrangianIntegrator_InitSolve, self );

    sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", Energy_SLE, False, NULL );
    if( sle ) {
        /* also set sle to run where required */
        EP_InsertClassHookAfter( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "SemiLagrangianIntegrator_InitSolve", 
            SystemLinearEquations_GetRunEPFunction(), sle );
        /* remember to disable the standard run at execute */
        SystemLinearEquations_SetRunDuringExecutePhase( sle, False );

        /* add the time step function */
        if( strcmp( sle->type, Energy_SLE_Type ) == 0 ) {
            EntryPoint_AppendClassHook( self->context->calcDtEP, "SemiLagrangianIntegrator_CalcAdvDiffDt", SemiLagrangianIntegrator_CalcAdvDiffDt, self->type, self );
        }
    }

    self->courant = Dictionary_GetDouble_WithDefault( self->context->dictionary, "courantFactor", 0.5 );

    self->isConstructed = True;
}

void _SemiLagrangianIntegrator_Build( void* slIntegrator, void* data ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i;

    if( self->velocityField ) Stg_Component_Build( self->velocityField, data, False );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Build( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Build( feVarStar , data, False );
    }
}

void _SemiLagrangianIntegrator_Initialise( void* slIntegrator, void* data ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i, el_i;
    unsigned 			patchSize 	= ( Mesh_GetDimSize( self->velocityField->feMesh ) == 2 ) ? 16 : 64;
    unsigned			nInc;

    self->dim = Mesh_GetDimSize( self->velocityField->feMesh );

    self->lElSize = Mesh_GetDomainSize( self->velocityField->feMesh, self->dim );

    if( self->velocityField ) Stg_Component_Initialise( self->velocityField, data, False );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Initialise( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Initialise( feVarStar , data, False );
    }

    self->ptsX = Memory_Alloc_2DArray_Unnamed( double, 4, 3 );
    self->ptsY = Memory_Alloc_2DArray_Unnamed( double, 4, 3 );
    self->ptsZ = Memory_Alloc_2DArray_Unnamed( double, 4, 3 );

    self->iArray = IArray_New();

    FeMesh_GetElementNodes( self->velocityField->feMesh, 0, self->iArray );
    nInc = IArray_GetSize( self->iArray );
    self->isQuad = ( nInc%3==0 ) ? True : False;

    if( !self->isQuad ) {
        self->elPatch = malloc( self->lElSize*sizeof(unsigned*) );
        for( el_i = 0; el_i < self->lElSize; el_i++ ) {
            self->elPatch[el_i] = malloc( patchSize*sizeof(unsigned) );
        }
        self->elPatchQuad = NULL;
        SemiLagrangianIntegrator_InitPatches( self );
    }
    else {
        unsigned nSections = ( self->dim == 2 ) ? 4 : 8;
        unsigned section_i;

        self->elPatch = NULL;
        self->elPatchQuad = malloc( self->lElSize*sizeof(unsigned**) );
        for( el_i = 0; el_i < self->lElSize; el_i++ ) {
            self->elPatchQuad[el_i] = malloc( nSections*sizeof(unsigned*) );
            for( section_i = 0; section_i < nSections; section_i++ ) {
                self->elPatchQuad[el_i][section_i] = malloc( patchSize*sizeof(unsigned) );
            }
        }
        SemiLagrangianIntegrator_InitPatches_Quad( self );
    }
}

void _SemiLagrangianIntegrator_Execute( void* slIntegrator, void* data ) {}

void _SemiLagrangianIntegrator_Destroy( void* slIntegrator, void* data ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i, el_i;

    if( !self->isQuad ) {
        for( el_i = 0; el_i < self->lElSize; el_i++ ) {
            free( self->elPatch[el_i] );
        }
        free( self->elPatch );
    }
    else {
        unsigned nSections = ( self->dim == 2 ) ? 4 : 8;
        unsigned section_i;

        for( el_i = 0; el_i < self->lElSize; el_i++ ) {
            for( section_i = 0; section_i < nSections; section_i++ ) {
                free( self->elPatchQuad[el_i][section_i] );
            }
            free( self->elPatchQuad[el_i] );
        }
        free( self->elPatchQuad );
    }

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Destroy( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Destroy( feVarStar , data, False );
    }
    if( self->pureAdvection ) Memory_Free( self->pureAdvection );

    Memory_Free( self->ptsX );
    Memory_Free( self->ptsY );
    Memory_Free( self->ptsZ );

    if( self->velocityField ) Stg_Component_Destroy( self->velocityField, data, False );

    Stg_Class_Delete( self->iArray );
}

void SemiLagrangianIntegrator_InitSolve( void* _self, void* _context ) {
    SemiLagrangianIntegrator*	self			= (SemiLagrangianIntegrator*) _self;
    unsigned			field_i, node_i;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    double            		phi[3];

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];

        /* generate the _* field */
        SemiLagrangianIntegrator_Solve( self, feVariable, feVarStar );

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
void SemiLagrangianIntegrator_IntegrateRK4( void* slIntegrator, double dt, double* origin, double* position ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeVariable*			velocityField	= self->velocityField;
    unsigned			ndims		= self->dim;
    unsigned			dim_i;
    double			min[3], max[3];
    double			k[4][3];
    double			coordPrime[3];
    unsigned*			periodic	= ((CartesianGenerator*)velocityField->feMesh->generator)->periodic;

    Mesh_GetGlobalCoordRange( velocityField->feMesh, min, max );

    SemiLagrangianIntegrator_CubicInterpolator( self, velocityField, origin, k[0] );
    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5 * dt * k[0][dim_i];
    }
    if( SemiLagrangianIntegrator_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SemiLagrangianIntegrator_CubicInterpolator( self, velocityField, coordPrime, k[1] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5 * dt * k[1][dim_i];
    }
    if( SemiLagrangianIntegrator_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SemiLagrangianIntegrator_CubicInterpolator( self, velocityField, coordPrime, k[2] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - dt * k[2][dim_i];
    }
    if( SemiLagrangianIntegrator_PeriodicUpdate( coordPrime, min, max, periodic, ndims ) );
    SemiLagrangianIntegrator_CubicInterpolator( self, velocityField, coordPrime, k[3] );

    for( dim_i = 0; dim_i < ndims; dim_i++ ) {
        position[dim_i] = origin[dim_i] -
            INV6 * dt * ( k[0][dim_i] + 2.0 * k[1][dim_i] + 2.0 * k[2][dim_i] + k[3][dim_i] );
    }
    if( SemiLagrangianIntegrator_PeriodicUpdate( position, min, max, periodic, ndims ) );
}

/* cubic Lagrangian interpoation in 1-D */
void SemiLagrangianIntegrator_InterpLagrange( double x, double* coords, double** values, unsigned numdofs, double* result ) {
    unsigned	node_i, dof_i;
    unsigned	otherIndices[3];
    unsigned	otherIndexCount, otherIndex_i;
    double	factor;

    for( dof_i = 0; dof_i < numdofs; dof_i++ )
        result[dof_i] = 0.0;

    for( node_i = 0; node_i < 4; node_i++ ) {
        otherIndexCount = 0;
        for( otherIndex_i = 0; otherIndex_i < 4; otherIndex_i++ )
            if( otherIndex_i != node_i )
                otherIndices[otherIndexCount++] = otherIndex_i;

        factor = 1.0;
        for( otherIndex_i = 0; otherIndex_i < 3; otherIndex_i++ )
            factor *= ( x - coords[otherIndices[otherIndex_i]] ) / ( coords[node_i] - coords[otherIndices[otherIndex_i]] );

        for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
            result[dof_i] += ( values[node_i][dof_i] * factor );
        }
    }
}

Bool SemiLagrangianIntegrator_PeriodicUpdate( double* pos, double* min, double* max, Bool* isPeriodic, unsigned ndim ) {
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

void SemiLagrangianIntegrator_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* pos, double* result ) {
    SemiLagrangianIntegrator*	self 			= (SemiLagrangianIntegrator*)slIntegrator;
    FeMesh*			feMesh			= feVariable->feMesh;
    int				dim			= self->dim;
    double			px[4], py[4], pz[4];
    unsigned			nDims			= Mesh_GetDimSize( feMesh );
    unsigned			numdofs			= feVariable->dofLayout->dofCounts[0];
    unsigned			elInd, index = 0;
    unsigned			x_i, y_i, z_i;
    unsigned*			patch;

    Mesh_SearchElements( feMesh, pos, &elInd );

    if( !self->isQuad ) {
        patch = self->elPatch[elInd];
    }
    else {
        int* inc;
        double* mid;

        FeMesh_GetElementNodes( feMesh, elInd, self->iArray );
        inc = IArray_GetPtr( self->iArray );        
        if( dim == 3 ) {
            mid = Mesh_GetVertex( feMesh, inc[13] );
            if( pos[0] < mid[0] && pos[1] < mid[1] && pos[2] < mid[2] )
                patch = self->elPatchQuad[elInd][0];
            else if( pos[0] > mid[0] && pos[1] < mid[1] && pos[2] < mid[2] )
                patch = self->elPatchQuad[elInd][1];
            else if( pos[0] < mid[0] && pos[1] > mid[1] && pos[2] < mid[2] )
                patch = self->elPatchQuad[elInd][2];
            else if( pos[0] > mid[0] && pos[1] > mid[1] && pos[2] < mid[2] )
                patch = self->elPatchQuad[elInd][3];
            else if( pos[0] < mid[0] && pos[1] < mid[1] && pos[2] > mid[2] )
                patch = self->elPatchQuad[elInd][4];
            else if( pos[0] > mid[0] && pos[1] < mid[1] && pos[2] > mid[2] )
                patch = self->elPatchQuad[elInd][5];
            else if( pos[0] < mid[0] && pos[1] > mid[1] && pos[2] > mid[2] )
                patch = self->elPatchQuad[elInd][6];
            else
                patch = self->elPatchQuad[elInd][7];
        }
        else {
            mid = Mesh_GetVertex( feMesh, inc[4] );
            if( pos[0] < mid[0] && pos[1] < mid[1] )
                patch = self->elPatchQuad[elInd][0];
            else if( pos[0] > mid[0] && pos[1] < mid[1] )
                patch = self->elPatchQuad[elInd][1];
            else if( pos[0] < mid[0] && pos[1] > mid[1] )
                patch = self->elPatchQuad[elInd][2];
            else
                patch = self->elPatchQuad[elInd][3];
        }
    }

    for( x_i = 0; x_i < 4; x_i++ )
        px[x_i] = Mesh_GetVertex( feMesh, patch[x_i] )[0];
    for( y_i = 0; y_i < 4; y_i++ )
        py[y_i] = Mesh_GetVertex( feMesh, patch[4*y_i] )[1];

    if( nDims == 3 ) {
        for( z_i = 0; z_i < 4; z_i++ )
            pz[z_i] = Mesh_GetVertex( feMesh, patch[16*z_i] )[2];

        for( z_i = 0; z_i < 4; z_i++ ) {
            for( y_i = 0; y_i < 4; y_i++ ) {
                for( x_i = 0; x_i < 4; x_i++ )
                    FeVariable_GetValueAtNode( feVariable, patch[index++], self->ptsX[x_i] );

                SemiLagrangianIntegrator_InterpLagrange( pos[0], px, self->ptsX, numdofs, self->ptsY[y_i] );
            }
            SemiLagrangianIntegrator_InterpLagrange( pos[1], py, self->ptsY, numdofs, self->ptsZ[z_i] );
        }
        SemiLagrangianIntegrator_InterpLagrange( pos[2], pz, self->ptsZ, numdofs, result );
    }
    else {
        for( y_i = 0; y_i < 4; y_i++ ) {
            for( x_i = 0; x_i < 4; x_i++ )
                FeVariable_GetValueAtNode( feVariable, patch[index++], self->ptsX[x_i] );

            SemiLagrangianIntegrator_InterpLagrange( pos[0], px, self->ptsX, numdofs, self->ptsY[y_i] );
        }
        SemiLagrangianIntegrator_InterpLagrange( pos[1], py, self->ptsY, numdofs, result );
    }
}

void SemiLagrangianIntegrator_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* varStarField ) {
    SemiLagrangianIntegrator*	self 		     = (SemiLagrangianIntegrator*)slIntegrator;
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

        SemiLagrangianIntegrator_IntegrateRK4( self, dt, coord, position );
        SemiLagrangianIntegrator_CubicInterpolator( self, variableField, position, var );
        FeVariable_SetValueAtNode( varStarField, node_I, var );
    }
    FeVariable_SyncShadowValues( varStarField );
}

Bool SemiLagrangianIntegrator_HasSide2D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes, int* sideNodes ) {
    int	nInc;

    Mesh_GetIncidence( feMesh, MT_VERTEX, elNodes[sideNodes[0]], MT_FACE, inc );
    nInc   = IArray_GetSize( inc );
    Mesh_GetIncidence( feMesh, MT_VERTEX, elNodes[sideNodes[1]], MT_FACE, inc );
    nInc  += IArray_GetSize( inc );
    if( nInc > 5 )
        return True;

    return False;
}

Bool SemiLagrangianIntegrator_HasLeft2D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (nNodes%3==0) ? 6 : 2;
    return SemiLagrangianIntegrator_HasSide2D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasRight2D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = (nNodes%3==0) ? 2 : 1;
    sideNodes[1] = (nNodes%3==0) ? 8 : 3;
    return SemiLagrangianIntegrator_HasSide2D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasBottom2D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (nNodes%3==0) ? 2 : 1;
    return SemiLagrangianIntegrator_HasSide2D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasTop2D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = (nNodes%3==0) ? 6 : 2;
    sideNodes[1] = (nNodes%3==0) ? 8 : 3;
    return SemiLagrangianIntegrator_HasSide2D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasSide3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes, int* sideNodes ) {
    int	nInc = 0, ii;

    for( ii = 0; ii < 4; ii++ ) {
        Mesh_GetIncidence( feMesh, MT_VERTEX, elNodes[sideNodes[ii]], Mesh_GetDimSize( feMesh ), inc );
        nInc += IArray_GetSize( inc );
    }
    if( nInc > 17 )
        return True;

    return False;
}

Bool SemiLagrangianIntegrator_HasLeft3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (nNodes%3==0) ?  6 : 2;
    sideNodes[2] = (nNodes%3==0) ? 18 : 4;
    sideNodes[3] = (nNodes%3==0) ? 24 : 6;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasRight3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = (nNodes%3==0) ?  2 : 1;
    sideNodes[1] = (nNodes%3==0) ?  8 : 3;
    sideNodes[2] = (nNodes%3==0) ? 20 : 5;
    sideNodes[3] = (nNodes%3==0) ? 26 : 7;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasBottom3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (nNodes%3==0) ?  2 : 1;
    sideNodes[2] = (nNodes%3==0) ? 18 : 4;
    sideNodes[3] = (nNodes%3==0) ? 20 : 5;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasTop3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = (nNodes%3==0) ?  6 : 2;
    sideNodes[1] = (nNodes%3==0) ?  8 : 3;
    sideNodes[2] = (nNodes%3==0) ? 24 : 6;
    sideNodes[3] = (nNodes%3==0) ? 26 : 7;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasFront3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = 0;
    sideNodes[1] = (nNodes%3==0) ? 2 : 1;
    sideNodes[2] = (nNodes%3==0) ? 6 : 2;
    sideNodes[3] = (nNodes%3==0) ? 8 : 3;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

Bool SemiLagrangianIntegrator_HasBack3D( FeMesh* feMesh, IArray* inc, int elInd, int* elNodes, int nNodes ) {
    int	sideNodes[4];

    sideNodes[0] = (nNodes%3==0) ? 18 : 4;
    sideNodes[1] = (nNodes%3==0) ? 20 : 5;
    sideNodes[2] = (nNodes%3==0) ? 24 : 6;
    sideNodes[3] = (nNodes%3==0) ? 26 : 7;

    return SemiLagrangianIntegrator_HasSide3D( feMesh, inc, elInd, elNodes, nNodes, sideNodes );
}

void SemiLagrangianIntegrator_InitPatches( void* slIntegrator ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    Grid*			grid		= ((CartesianGenerator*)feMesh->generator)->vertGrid;
    int				nInc, *inc;
    unsigned			lNode, gNode, el_i;
    unsigned			index, IJK[3], IJK_test[3];
    IArray*			iArray1		= IArray_New();
    IArray*			iArray2		= IArray_New();

    for( el_i = 0; el_i < Mesh_GetDomainSize( feMesh, self->dim ); el_i++ ) {
        FeMesh_GetElementNodes( feMesh, el_i, iArray1 );
        inc  = IArray_GetPtr( iArray1 );
        nInc = IArray_GetSize( iArray1 );

        gNode = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
        Grid_Lift( grid, gNode, IJK );

        if( self->dim == 3 ) {
            if( !SemiLagrangianIntegrator_HasRight3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( !SemiLagrangianIntegrator_HasBack3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;

            if( SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
        }
        else {
            if( !SemiLagrangianIntegrator_HasRight2D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop2D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;

            if( SemiLagrangianIntegrator_HasLeft2D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( SemiLagrangianIntegrator_HasBottom2D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        }

        FeMesh_NodeGlobalToDomain( feMesh, Grid_Project( grid, IJK ), &lNode );

        index = 0;
        /* interpolate using Lagrange polynomial shape functions */
        if( self->dim == 3 ) {
            for( IJK_test[2] = IJK[2]; IJK_test[2] < IJK[2] + 4; IJK_test[2]++ ) {
                for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
                    for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
                        gNode = Grid_Project( grid, IJK_test );
                        if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode, &lNode ) ) abort();
                        else
                            self->elPatch[el_i][index++] = lNode;
                    }
                }
            }
        }
        else {
            for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
                for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
                    gNode = Grid_Project( grid, IJK_test );
                    if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode, &lNode ) ) abort();
                    else
                        self->elPatch[el_i][index++] = lNode;
                }
            }
        }
    }

    Stg_Class_Delete( iArray1 );
    Stg_Class_Delete( iArray2 );
}

void SemiLagrangianIntegrator_InitPatches_Quad( void* slIntegrator ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    Grid*			grid		= ((CartesianGenerator*)feMesh->generator)->vertGrid;
    IArray*			iArray1		= IArray_New();
    IArray*			iArray2		= IArray_New();
    unsigned			el_i, lNode, gNode;
    unsigned			IJK[3], IJK_test[3], index;
    int				*inc, nInc;
    int				gOrig[8], orig_i;

    for( el_i = 0; el_i < Mesh_GetDomainSize( feMesh, self->dim ); el_i++ ) {
        FeMesh_GetElementNodes( feMesh, el_i, iArray1 );
        inc  = IArray_GetPtr( iArray1 );
        nInc = IArray_GetSize( iArray1 );

        if( self->dim == 3 ) {
            gNode = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
            /* left-bottom-front */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[0] = Grid_Project( grid, IJK );
            /* right-bottom-front */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[1] = Grid_Project( grid, IJK );
            /* left-top-front */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[2] = Grid_Project( grid, IJK );
            /* right-top-front */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[3] = Grid_Project( grid, IJK );
            /* left-bottom-back */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[4] = Grid_Project( grid, IJK );
            /* right-bottom-back */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[5] = Grid_Project( grid, IJK );
            /* left-top-back */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[6] = Grid_Project( grid, IJK );
            /* right-top-back */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[2]--;
            gOrig[7] = Grid_Project( grid, IJK );

            for( orig_i = 0; orig_i < 8; orig_i++ ) {
                Grid_Lift( grid, gOrig[orig_i], IJK );
                FeMesh_NodeGlobalToDomain( feMesh, gOrig[orig_i], &lNode );
                index = 0;
                for( IJK_test[2] = IJK[2]; IJK_test[2] < IJK[2] + 4; IJK_test[2]++ ) {
                    for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
                        for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
                            gNode = Grid_Project( grid, IJK_test );
                            if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode, &lNode ) ) abort();
                            else
                                self->elPatchQuad[el_i][orig_i][index++] = lNode;
                        }
                    }
                }
            }
        }
        else {
            gNode = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
            /* left-bottom-front */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft2D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom2D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            gOrig[0] = Grid_Project( grid, IJK );
            /* right-bottom-front */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight2D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if(  SemiLagrangianIntegrator_HasBottom2D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            gOrig[1] = Grid_Project( grid, IJK );
            /* left-top-front */
            Grid_Lift( grid, gNode, IJK );
            if(  SemiLagrangianIntegrator_HasLeft2D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop2D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            gOrig[2] = Grid_Project( grid, IJK );
            /* right-top-front */
            Grid_Lift( grid, gNode, IJK );
            if( !SemiLagrangianIntegrator_HasRight2D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
            if( !SemiLagrangianIntegrator_HasTop2D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
            gOrig[3] = Grid_Project( grid, IJK );

            for( orig_i = 0; orig_i < 4; orig_i++ ) {
                Grid_Lift( grid, gOrig[orig_i], IJK );
                FeMesh_NodeGlobalToDomain( feMesh, gOrig[orig_i], &lNode );
                index = 0;
                for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
                    for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
                        gNode = Grid_Project( grid, IJK_test );
                        if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode, &lNode ) ) abort();
                        else
                            self->elPatchQuad[el_i][orig_i][index++] = lNode;
                    }
                }
            }
        }
    }

    Stg_Class_Delete( iArray1 );
    Stg_Class_Delete( iArray2 );
}
#if 0
void SemiLagrangianIntegrator_InitPatches( void* slIntegrator ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    double			min[3], max[3], delta[3], pos[3];
    unsigned			dim		= Mesh_GetDimSize( feMesh );
    unsigned			el_i, *inc, gNode_i, lNode_i, nNodes[3], index;
    int				x_0, y_0, z_0, x_i, y_i, z_i;

    SemiLagrangianIntegrator_GetDeltaConst( self->velocityField, delta, nNodes );

    for( el_i = 0; el_i < Mesh_GetLocalSize( feMesh, dim ); el_i++ ) {
        FeMesh_GetElementNodes( feMesh, el_i, self->velocityField->inc );
        inc  = IArray_GetPtr( self->velocityField->inc );

        Mesh_GetLocalCoordRange( feMesh, min, max );
        gNode_i = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );

        x_0 = (int) gNode_i % nNodes[0];
        y_0 = ((int) gNode_i / nNodes[0])%nNodes[1];
        if( dim == 3 )
            z_0 = (int) gNode_i / ( nNodes[0]*nNodes[1] );

        pos[0] = Mesh_GetVertex( feMesh, inc[0] )[0] + 1.0e-4;
        pos[1] = Mesh_GetVertex( feMesh, inc[0] )[1] + 1.0e-4;
        if( dim == 3 )
            pos[2] = Mesh_GetVertex( feMesh, inc[0] )[2] + 1.0e-4;

        if( pos[0] >  min[0] + delta[0] ) x_0--;
        if( pos[0] >= max[0] - delta[0] ) x_0--;

        if( pos[1] >  min[1] + delta[1] ) y_0--;
        if( pos[1] >= max[1] - delta[1] ) y_0--;

        if( dim == 3 ) {
            if( pos[2] >  min[2] + delta[2] ) z_0--;
            if( pos[2] >= max[2] - delta[2] ) z_0--;
        }

        index = 0;
        if( dim == 2 ) {
            for( y_i = 0; y_i < 4; y_i++ )
                for( x_i = 0; x_i < 4; x_i++ ) {
                    gNode_i = x_0 + x_i + ( y_0 + y_i ) * nNodes[0];
                    if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_i, &lNode_i ) ) abort();
                    else
                        self->elPatch[el_i][index++] = lNode_i;
                }
        }
        else {
            for( z_i = 0; z_i < 4; z_i++ )
                for( y_i = 0; y_i < 4; y_i++ )
                    for( x_i = 0; x_i < 4; x_i++ ) {
                        gNode_i = x_0 + x_i + ( y_0 + y_i ) * nNodes[0] + ( z_0 + z_i ) * nNodes[0] * nNodes[1];
                        if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_i, &lNode_i ) ) abort();
                        else
                            self->elPatch[el_i][index++] = lNode_i;
                    }
        }
    }
}
#endif

double SemiLagrangianIntegrator_CalcAdvDiffDt( void* slIntegrator, FiniteElementContext* context ) {
    SemiLagrangianIntegrator*	self 		= (SemiLagrangianIntegrator*)slIntegrator;
    double			lAdv, lDif, gAdv, gDif, dt, dx[3], dxMin, vMag;
    double 			manualDt        = Dictionary_GetDouble_WithDefault( context->dictionary, "manualTimeStep", 0.0 );

    if( manualDt > 1.0e-6 ) {
        return manualDt;
    }

    FeVariable_GetMinimumSeparation( self->velocityField, &dxMin, dx );
    vMag = FieldVariable_GetMaxGlobalFieldMagnitude( self->velocityField );
    lAdv = self->courant*dxMin/vMag;
    lDif = self->courant*dxMin*dxMin;

    MPI_Allreduce( &lAdv, &gAdv, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &lDif, &gDif, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

    dt = ( gAdv < gDif ) ? gAdv : gDif;

    return dt;
}
