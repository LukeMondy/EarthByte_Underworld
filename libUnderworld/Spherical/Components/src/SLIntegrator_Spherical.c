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

/** Textual name of this class */
const Type SLIntegrator_Spherical_Type = "SLIntegrator_Spherical";

SLIntegrator_Spherical* _SLIntegrator_Spherical_New( SLINTEGRATOR_SPHERICAL_DEFARGS ) {
    SLIntegrator_Spherical*		self;

    /* Allocate memory */
    assert( _sizeOfSelf >= sizeof(SLIntegrator_Spherical) );
    /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
    /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
       and so should be set to ZERO in any children of this class. */
    nameAllocationType = NON_GLOBAL;

    self = (SLIntegrator_Spherical*) _Stg_Component_New( STG_COMPONENT_PASSARGS );

    /* General info */
    self->variableList = Stg_ObjectList_New();
    self->varStarList  = Stg_ObjectList_New();

    return self;
}

void* _SLIntegrator_Spherical_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
    SLIntegrator_Spherical*	self 	= (SLIntegrator_Spherical*)slIntegrator;
    SLIntegrator_Spherical*	newSLIntegrator_Spherical;
    PtrMap*			map 	= ptrMap;
    Bool			ownMap 	= False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSLIntegrator_Spherical = _Stg_Component_Copy( self, dest, deep, nameExt, map );

    if( deep ) {
        if( (newSLIntegrator_Spherical->velocityField = PtrMap_Find( map, self->velocityField )) == NULL ) {
            newSLIntegrator_Spherical->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
            PtrMap_Append( map, self->velocityField, newSLIntegrator_Spherical->velocityField );
        }
    }
    else {
        newSLIntegrator_Spherical->velocityField = Stg_Class_Copy( self->velocityField, NULL, deep, nameExt, map );
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSLIntegrator_Spherical;
}


void _SLIntegrator_Spherical_Delete( void* slIntegrator ) {
    SLIntegrator_Spherical*		self = (SLIntegrator_Spherical*)slIntegrator;
    Stg_Class_Delete( self->variableList );
    Stg_Class_Delete( self->varStarList );
}

void _SLIntegrator_Spherical_Print( void* slIntegrator, Stream* stream ) {
    SLIntegrator_Spherical*		self = (SLIntegrator_Spherical*)slIntegrator;

    _Stg_Component_Print( self, stream );

    Journal_PrintPointer( stream, self->velocityField );
}

void* _SLIntegrator_Spherical_DefaultNew( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(SLIntegrator_Spherical);
    Type                                                      type = SLIntegrator_Spherical_Type;
    Stg_Class_DeleteFunction*                              _delete = _SLIntegrator_Spherical_Delete;
    Stg_Class_PrintFunction*                                _print = _SLIntegrator_Spherical_Print;
    Stg_Class_CopyFunction*                                  _copy = _SLIntegrator_Spherical_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SLIntegrator_Spherical_DefaultNew;
    Stg_Component_ConstructFunction*                    _construct = _SLIntegrator_Spherical_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = _SLIntegrator_Spherical_Build;
    Stg_Component_InitialiseFunction*                  _initialise = _SLIntegrator_Spherical_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = _SLIntegrator_Spherical_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _SLIntegrator_Spherical_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL; /* default value NON_GLOBAL */

    return (void*)_SLIntegrator_Spherical_New( SLINTEGRATOR_SPHERICAL_PASSARGS );
}

void _SLIntegrator_Spherical_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data ) {
    SLIntegrator_Spherical*	self 		= (SLIntegrator_Spherical*)slIntegrator;
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

    EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), SLIntegrator_Spherical_InitSolve, self );

    sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"SLE", Energy_SLE, False, NULL );
    if( sle ) {
        /* also set sle to run where required */
        EP_InsertClassHookAfter( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ), "SLIntegrator_Spherical_InitSolve", 
            SystemLinearEquations_GetRunEPFunction(), sle );
        /* remember to disable the standard run at execute */
        SystemLinearEquations_SetRunDuringExecutePhase( sle, False );

        /* add the time step function */
        if( strcmp( sle->type, Energy_SLE_Type ) == 0 ) {
            EntryPoint_AppendClassHook( self->context->calcDtEP, "SLIntegrator_Spherical_CalcAdvDiffDt", SLIntegrator_Spherical_CalcAdvDiffDt, SLIntegrator_Spherical_Type, self );
        }
    }

    self->courant = Dictionary_GetDouble_WithDefault( self->context->dictionary, "courantFactor", 0.5 );

    self->isConstructed = True;
}

void _SLIntegrator_Spherical_Build( void* slIntegrator, void* data ) {
    SLIntegrator_Spherical*	self 		= (SLIntegrator_Spherical*)slIntegrator;
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

void _SLIntegrator_Spherical_Initialise( void* slIntegrator, void* data ) {
    SLIntegrator_Spherical*	self 		= (SLIntegrator_Spherical*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    unsigned			field_i, el_i;

    if( self->velocityField ) Stg_Component_Initialise( self->velocityField, data, False );

    self->dElSize = Mesh_GetDomainSize( self->velocityField->feMesh, MT_VOLUME );

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];
        if( feVariable ) Stg_Component_Initialise( feVariable, data, False );
        if( feVarStar  ) Stg_Component_Initialise( feVarStar , data, False );
    }

    self->abcissa = malloc(4*sizeof(double));
    self->Ni      = malloc(64*sizeof(double));
    self->GNix    = malloc(MT_VOLUME*sizeof(double*));
    self->GNix[0] = malloc(64*sizeof(double));
    self->GNix[1] = malloc(64*sizeof(double));
    self->GNix[2] = malloc(64*sizeof(double));

    /* sanity checks */
    if( Mesh_GetDimSize( feMesh ) != MT_VOLUME ) {
        printf( "ERROR: component %s requires mesh of dimension 3\n", self->type );
        abort();
    }
    if( strcmp( feMesh->algorithms->type, "Mesh_SphericalAlgorithms" ) != 0 ) {
        printf( "ERROR: component %s required mesh algorithms type Mesh_SphericalAlgoritms, type is %s.\n", self->type, feMesh->algorithms->type );
        abort();
    }

    FeMesh_GetElementNodes( feMesh, 0, self->inc );
    self->elQuad = ( IArray_GetSize( self->inc )%3==0 ) ? True : False;

    if( !self->elQuad ) {
        self->elPatch = malloc(Mesh_GetDomainSize( feMesh, MT_VOLUME )*sizeof(unsigned*));
        for( el_i = 0; el_i < Mesh_GetDomainSize( feMesh, MT_VOLUME ); el_i++ ) {
            self->elPatch[el_i] = malloc( 64*sizeof(unsigned) );
        }
        self->elPatchQuad = NULL;
        SLIntegrator_Spherical_InitPatches( self );
    }
    else {
        int section_i;

        self->elPatch = NULL;
        self->elPatchQuad = malloc(Mesh_GetDomainSize( self->velocityField->feMesh, MT_VOLUME )*sizeof(unsigned**));
        for( el_i = 0; el_i < Mesh_GetDomainSize( self->velocityField->feMesh, MT_VOLUME ); el_i++ ) {
            self->elPatchQuad[el_i] = malloc( 8*sizeof(unsigned*) );
            for( section_i = 0; section_i < 4; section_i++ ) {
                self->elPatchQuad[el_i][section_i] = malloc( 16*sizeof(unsigned) );
            }
        }
        SLIntegrator_Spherical_InitPatches_Quad( self );
    }

    self->abcissa[0] = -1.0;
    self->abcissa[0] = -1.0;
    self->abcissa[1] = -0.33333333333333333;
    self->abcissa[2] = -1.0*self->abcissa[1];
    self->abcissa[3] = -1.0*self->abcissa[0];
}

void _SLIntegrator_Spherical_Execute( void* slIntegrator, void* data ) {}

void _SLIntegrator_Spherical_Destroy( void* slIntegrator, void* data ) {
    SLIntegrator_Spherical*	self 		= (SLIntegrator_Spherical*)slIntegrator;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    unsigned			field_i, el_i;

    if( !self->elQuad ) {
        for( el_i = 0; el_i < self->dElSize; el_i++ ) {
            free( self->elPatch[el_i] );
        }
        free( self->elPatch );
    }
    else {
        int section_i;

        for( el_i = 0; el_i < self->dElSize; el_i++ ) {
            for( section_i = 0; section_i < 8; section_i++ ) {
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

    free(self->abcissa);
    free(self->Ni);
    free(self->GNix[0]);
    free(self->GNix[1]);
    free(self->GNix[2]);
    free(self->GNix);

    Stg_Class_Delete( self->inc );

    if( self->velocityField ) Stg_Component_Destroy( self->velocityField, data, False );
}

void SLIntegrator_Spherical_InitSolve( void* _self, void* _context ) {
    SLIntegrator_Spherical*	self			= (SLIntegrator_Spherical*) _self;
    unsigned			field_i, node_i;
    FeVariable*			feVariable;
    FeVariable*			feVarStar;
    double            		phi[3];

    for( field_i = 0; field_i < self->variableList->count; field_i++ ) {
        feVariable = (FeVariable*) self->variableList->data[field_i];
        feVarStar  = (FeVariable*) self->varStarList->data[field_i];

        /* generate the _* field */
        SLIntegrator_Spherical_Solve( self, feVariable, feVarStar );

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
void SLIntegrator_Spherical_IntegrateRK4( void* slIntegrator, FeVariable* velocityField, double dt, double* origin, double* position ) {
    SLIntegrator_Spherical*	self 	     	= (SLIntegrator_Spherical*)slIntegrator;
    unsigned			dim_i;
    double			k[4][3];
    double			coordPrime[3];

    SLIntegrator_Spherical_CubicInterpolator( self, velocityField, origin, k[0] );
    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[0][dim_i];
    }
    SLIntegrator_Spherical_BoundaryUpdate3D( velocityField->feMesh, self->inc, coordPrime );
    SLIntegrator_Spherical_CubicInterpolator( self, velocityField, coordPrime, k[1] );

    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[1][dim_i];
    }
    SLIntegrator_Spherical_BoundaryUpdate3D( velocityField->feMesh, self->inc, coordPrime );
    SLIntegrator_Spherical_CubicInterpolator( self, velocityField, coordPrime, k[2] );

    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        coordPrime[dim_i] = origin[dim_i] - dt*k[2][dim_i];
    }
    SLIntegrator_Spherical_BoundaryUpdate3D( velocityField->feMesh, self->inc, coordPrime );
    SLIntegrator_Spherical_CubicInterpolator( self, velocityField, coordPrime, k[3] );

    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        position[dim_i] = origin[dim_i] - INV6*dt*( k[0][dim_i] + 2.0*k[1][dim_i] + 2.0*k[2][dim_i] + k[3][dim_i] );
    }
    SLIntegrator_Spherical_BoundaryUpdate3D( velocityField->feMesh, self->inc, position );
}

void Spherical_XYZ2regionalSphere( double *xyz, double* rs ) {
    rs[0] = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
    rs[1] = atan2( xyz[0], xyz[2] );
    rs[2] = atan2( xyz[1], xyz[2] );
}

void Spherical_RegionalSphere2XYZ( double *rs, double* xyz ) {
    double r 		= rs[0];
    double theta 	= rs[1];
    double phi 		= rs[2];
    double X 		= tan(theta);
    double Y		= tan(phi);
    double d 		= sqrt(1 + X*X + Y*Y);

    xyz[0] = r/d * X;
    xyz[1] = r/d * Y;
    xyz[2] = r/d * 1;
}

void Spherical_VectorXYZ2regionalSphere( double* v, double* x, double* v2 ) {
    double	rInv	= 1.0/sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
    double	a2	= x[0]*x[0]/x[2]/x[2];
    double	b2	= x[1]*x[1]/x[2]/x[2];
    double	dr[3], dEta[3], dZeta[3];
    unsigned	i;

    dr[0]    = x[0]*rInv;
    dr[1]    = x[1]*rInv;
    dr[2]    = x[2]*rInv;
    dEta[0]  = 1.0/x[2]/(1.0 + a2);
    dEta[1]  = 0.0;
    dEta[2]  = -x[0]/x[2]/x[2]/(1.0 + a2);
    dZeta[0] = 0.0;
    dZeta[1] = 1.0/x[2]/(1.0 + b2);
    dZeta[2] = -x[1]/x[2]/x[2]/(1.0 + b2);

    v2[0] = v2[1] = v2[2] = 0.0;

    for( i = 0; i < 3; i++ ) {
        v2[0] += dr[i]*v[i];
        v2[1] += dEta[i]*v[i];
        v2[2] += dZeta[i]*v[i];
    }
}

#define SL_EPS 1.0e-04
void SLIntegrator_Spherical_BoundaryUpdate3D( FeMesh* feMesh, IArray* iArray, double* pos ) {
    unsigned*   periodic        = ((CartesianGenerator*)feMesh->generator)->periodic;
    double      rs[3], min[3], max[3];
    Grid**      grid            = (Grid**) Mesh_GetExtension( feMesh, Grid*, feMesh->elGridId );
    unsigned*   sizes           = Grid_GetSizes( *grid );
    unsigned	dim_i;

    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        min[dim_i] = ((RSGenerator*)feMesh->generator)->crdMin[dim_i];
        max[dim_i] = ((RSGenerator*)feMesh->generator)->crdMax[dim_i];
    }
    for( dim_i = 1; dim_i < MT_VOLUME; dim_i++ ) {
        min[dim_i] *= M_PI/180;
        max[dim_i] *= M_PI/180;
    }

    Spherical_XYZ2regionalSphere( pos, rs );

    /*if( !periodic[0] ) {
        double dr = (max[0]-min[0])/sizes[0];
        double testRS[3], testXYZ[3], v1[3], v2[3], norm[3], nMag, hyp[3], hNorm, lambda, theta;
        unsigned nInc, *inc, ind0, ind1, ind2;
        if( rs[0] > max[0] - SL_EPS ) {
            testRS[0] = rs[0] - 0.5*dr;
            for( dim_i = 1; dim_i < MT_VOLUME; dim_i++ ) {
                testRS[dim_i] = rs[dim_i];
                if( testRS[dim_i] > max[dim_i] ) testRS[dim_i] = max[dim_i];
                if( testRS[dim_i] < min[dim_i] ) testRS[dim_i] = min[dim_i];
            }
            Spherical_RegionalSphere2XYZ( testRS, testXYZ );
            Mesh_SearchElements( feMesh, testXYZ, &elInd );
            FeMesh_GetElementNodes( feMesh, elInd, iArray );
            nInc = IArray_GetSize( iArray );
            inc  = IArray_GetPtr( iArray );
            if( nInc%3 == 0 ) {
                ind0 = 2; ind1 = 8; ind2 = 20;
            } else {
                ind0 = 1; ind1 = 3; ind2 = 5;
            }
            for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
                v1[dim_i] = Mesh_GetVertex( feMesh, inc[ind1] )[dim_i] - Mesh_GetVertex( feMesh, inc[ind0] )[dim_i];
                v2[dim_i] = Mesh_GetVertex( feMesh, inc[ind2] )[dim_i] - Mesh_GetVertex( feMesh, inc[ind0] )[dim_i];
            }
            norm[0] = v1[1]*v2[2] - v1[2]*v2[1];
            norm[1] = v1[2]*v2[0] - v1[0]*v2[2];
            norm[2] = v1[0]*v2[1] - v1[1]*v2[0];
            nMag = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
            norm[0] /= nMag; norm[1] /= nMag; norm[2] /= nMag;
            for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
                hyp[dim_i] = pos[dim_i] - Mesh_GetVertex( feMesh, inc[ind0] )[dim_i];
            }
            hNorm  = sqrt( hyp[0]*hyp[0] + hyp[1]*hyp[1] + hyp[2]*hyp[2] );
            theta  = 0.25*M_PI - acos( ( norm[0]*hyp[0] + norm[1]*hyp[1] + norm[2]*hyp[2] )/hNorm );
            lambda = fabs( hNorm*sin( theta ) );
            for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
                pos[dim_i] -= lambda*norm[dim_i];
            }
            Spherical_XYZ2regionalSphere( pos, rs );

        } else if( rs[0] < min[0] ) {
            rs[0] = min[0];    
        }
    }*/

    for( dim_i = 0; dim_i < MT_VOLUME; dim_i++ ) {
        if( rs[dim_i] < min[dim_i] ) {
            rs[dim_i] = (periodic[dim_i]) ? max[dim_i] - min[dim_i] + rs[dim_i] : min[dim_i];
        }
        if( rs[dim_i] > max[dim_i] ) {
            rs[dim_i] = (periodic[dim_i]) ? min[dim_i] - max[dim_i] + rs[dim_i] : max[dim_i];
        }
    }
    Spherical_RegionalSphere2XYZ( rs, pos );
}

void SLIntegrator_Spherical_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* position, double* result ) {
    SLIntegrator_Spherical*	self 			= (SLIntegrator_Spherical*)slIntegrator;
    FeMesh*			feMesh			= feVariable->feMesh;
    unsigned			elInd;
    unsigned			numdofs			= feVariable->dofLayout->dofCounts[0];
    double			lCoord[3], phi_i[3];
    unsigned			node_i, dof_i;
    Grid**      		grid                    = (Grid**)Mesh_GetExtension( feMesh, Grid*, feMesh->vertGridId );
    unsigned*			patch;

    Mesh_SearchElements( feMesh, position, &elInd );

    if( !self->elQuad ) {
        patch = self->elPatch[elInd];
    }
    else {
        int* inc;
        double mid[3], rez[3];

        FeMesh_GetElementNodes( feMesh, elInd, self->inc );
        inc = IArray_GetPtr( self->inc );
        
        Spherical_XYZ2regionalSphere( Mesh_GetVertex( feMesh, inc[13] ), mid );
        Spherical_XYZ2regionalSphere( position, rez );

        if( rez[0] < mid[0] && rez[1] < mid[1] && rez[2] < mid[2] )
            patch = self->elPatchQuad[elInd][0];
        else if( rez[0] > mid[0] && rez[1] < mid[1] && rez[2] < mid[2] )
            patch = self->elPatchQuad[elInd][1];
        else if( rez[0] < mid[0] && rez[1] > mid[1] && rez[2] < mid[2] )
            patch = self->elPatchQuad[elInd][2];
        else if( rez[0] > mid[0] && rez[1] > mid[1] && rez[2] < mid[2] )
            patch = self->elPatchQuad[elInd][3];
        else if( rez[0] < mid[0] && rez[1] < mid[1] && rez[2] > mid[2] )
            patch = self->elPatchQuad[elInd][4];
        else if( rez[0] > mid[0] && rez[1] < mid[1] && rez[2] > mid[2] )
            patch = self->elPatchQuad[elInd][5];
        else if( rez[0] < mid[0] && rez[1] > mid[1] && rez[2] > mid[2] )
            patch = self->elPatchQuad[elInd][6];
        else
            patch = self->elPatchQuad[elInd][7];
    }

    SLIntegrator_Spherical_GlobalToLocal( self, feMesh, patch, position, lCoord );
    SLIntegrator_Spherical_ShapeFuncs3D( self, lCoord, self->Ni );
    for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
        result[dof_i] = 0.0;
    }
    for( node_i = 0; node_i < 64; node_i++ ) {
        FeVariable_GetValueAtNode( feVariable, patch[node_i], phi_i );
        for( dof_i = 0; dof_i < numdofs; dof_i++ ) {
            result[dof_i] += phi_i[dof_i]*self->Ni[node_i];
        }
    }
}

void SLIntegrator_Spherical_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* varStarField ) {
    SLIntegrator_Spherical*	self 		     = (SLIntegrator_Spherical*)slIntegrator;
    FiniteElementContext*	context		     = self->context;
    unsigned			node_i, gNode_i;
    FeMesh*			feMesh		     = variableField->feMesh;
    unsigned			meshSize	     = Mesh_GetLocalSize( feMesh, MT_VERTEX );
    FeVariable*			velocityField	     = self->velocityField;
    double			dt		     = AbstractContext_Dt( context );
    double			position[3];
    double			var[3];
    double*			coord;
    Grid**      		grid                 = (Grid**) Mesh_GetExtension( feMesh, Grid*, feMesh->vertGridId );
    unsigned*   		sizes                = Grid_GetSizes( *grid );

    FeVariable_SyncShadowValues( velocityField );
    FeVariable_SyncShadowValues( variableField );

    /* assume that the variable mesh is the same as the velocity mesh */
    for( node_i = 0; node_i < meshSize; node_i++ ) {
        gNode_i = Mesh_DomainToGlobal( feMesh, MT_VERTEX, node_i );
        /* skip if node is an outer wall with a boundary condition to avoid overshooting of characteristics 
           to departure points outside circle. TODO: only consistent for scalar fields at present */
        if( FeVariable_IsBC( variableField, node_i, 0 ) && gNode_i%sizes[0] == sizes[0] - 1 /*outer (right) wall*/ )
            continue;

	coord = Mesh_GetVertex( feMesh, node_i );

        SLIntegrator_Spherical_IntegrateRK4( self, velocityField, dt, coord, position );
        SLIntegrator_Spherical_CubicInterpolator( self, variableField, position, var );
        FeVariable_SetValueAtNode( varStarField, node_i, var );
    }
    FeVariable_SyncShadowValues( varStarField );
}

double SLIntegrator_Spherical_Lagrange( double* xi, double x, int j ) {
    int i;
    double l = 1.0;

    for( i = 0; i < 4; i++ ) {
        if( i == j )
            continue;

        l *= (x - xi[i])/(xi[j]-xi[i]);
    }
    return l;
}

double SLIntegrator_Spherical_LagrangeDeriv( double* xi, double x, int j ) {
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

void SLIntegrator_Spherical_ShapeFuncs3D( void* slIntegrator, double* xi, double* const Ni ) {
    SLIntegrator_Spherical* self = (SLIntegrator_Spherical*) slIntegrator;
    int                     x_j, y_j, z_j, pt_j = 0;
    double                  Nx, Ny, Nz;

    for( z_j = 0; z_j < 4; z_j++ ) {
        Nz = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[2], z_j );

        for( y_j = 0; y_j < 4; y_j++ ) {
            Ny = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[1], y_j );
            
            for( x_j = 0; x_j < 4; x_j++ ) {
                Nx = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[0], x_j );

                Ni[pt_j++] = Nx*Ny*Nz;
            }
        }
    }
}

void SLIntegrator_Spherical_ShapeFuncDerivs3D( void* slIntegrator, double* xi, double** const GNix ) {
    SLIntegrator_Spherical* self = (SLIntegrator_Spherical*) slIntegrator;
    double                  Nx, Ny, Nz, GNx, GNy, GNz;
    int                     x_j, y_j, z_j, pt_j = 0;

    for( z_j = 0; z_j < 4; z_j++ ) {
        Nz  = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[2], z_j );
        GNz = SLIntegrator_Spherical_LagrangeDeriv( self->abcissa, xi[2], z_j );

        for( y_j = 0; y_j < 4; y_j++ ) {
            Ny  = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[1], y_j );
            GNy = SLIntegrator_Spherical_LagrangeDeriv( self->abcissa, xi[1], y_j );
            
            for( x_j = 0; x_j < 4; x_j++ ) {
                Nx  = SLIntegrator_Spherical_Lagrange( self->abcissa, xi[0], x_j );
                GNx = SLIntegrator_Spherical_LagrangeDeriv( self->abcissa, xi[0], x_j );

                GNix[0][pt_j] = GNx*Ny*Nz;
                GNix[1][pt_j] = Nx*GNy*Nz;
                GNix[2][pt_j] = Nx*Ny*GNz;
                pt_j++;
            }
        }
    }
}

void SLIntegrator_Spherical_GlobalToLocal( void* slIntegrator, void* _mesh, unsigned* nodeInds, const double* gCoord, double* lCoord ) {
    SLIntegrator_Spherical* 	self 		= (SLIntegrator_Spherical*)slIntegrator;
    Mesh*			mesh 		= (Mesh*)_mesh;
    TensorArray         	jacobiMatrix;
    double              	tolerance       = 0.0001;
    double              	maxResidual;
    Iteration_Index     	maxIterations   = 100;
    Iteration_Index     	iteration_I;
    Node_Index          	node_I;
    Node_Index          	nodeCount       = 64;
    XYZ                 	rightHandSide;
    XYZ                 	xiIncrement;
    double*       	    	nodeCoord;
    double*		    	Ni		= self->Ni;
    double**	    		GNix		= self->GNix;

    /* Initial guess for element local coordinate is in the centre of the element - ( 0.0, 0.0, 0.0 ) */
    memset( lCoord, 0, MT_VOLUME*sizeof(double) );

    /* Do Newton-Raphson Iteration */
    for ( iteration_I = 0 ; iteration_I < maxIterations ; iteration_I++ ) {
        /* Initialise Values */
        TensorArray_Zero( jacobiMatrix );
        memset( rightHandSide, 0, sizeof( XYZ ) );

        /* Evaluate shape functions for rhs */
        SLIntegrator_Spherical_ShapeFuncs3D( self, lCoord, Ni );
        SLIntegrator_Spherical_ShapeFuncDerivs3D( self, lCoord, GNix );

        for( node_I = 0 ; node_I < nodeCount ; node_I++ ) {
            nodeCoord = Mesh_GetVertex( mesh, nodeInds[node_I] );

            /* Form jacobi matrix */
            jacobiMatrix[ MAP_3D_TENSOR( 0, 0 ) ] += GNix[0][node_I] * nodeCoord[ I_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 0, 1 ) ] += GNix[1][node_I] * nodeCoord[ I_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 0, 2 ) ] += GNix[2][node_I] * nodeCoord[ I_AXIS ];

            jacobiMatrix[ MAP_3D_TENSOR( 1, 0 ) ] += GNix[0][node_I] * nodeCoord[ J_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 1, 1 ) ] += GNix[1][node_I] * nodeCoord[ J_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 1, 2 ) ] += GNix[2][node_I] * nodeCoord[ J_AXIS ];

            jacobiMatrix[ MAP_3D_TENSOR( 2, 0 ) ] += GNix[0][node_I] * nodeCoord[ K_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 2, 1 ) ] += GNix[1][node_I] * nodeCoord[ K_AXIS ];
            jacobiMatrix[ MAP_3D_TENSOR( 2, 2 ) ] += GNix[2][node_I] * nodeCoord[ K_AXIS ];

            /* Form right hand side */
            rightHandSide[ I_AXIS ] -= Ni[node_I] * nodeCoord[ I_AXIS ];
            rightHandSide[ J_AXIS ] -= Ni[node_I] * nodeCoord[ J_AXIS ];
            rightHandSide[ K_AXIS ] -= Ni[node_I] * nodeCoord[ K_AXIS ];
        }

        /* Finish building right hand side */
        rightHandSide[ I_AXIS ] += gCoord[ I_AXIS ];
        rightHandSide[ J_AXIS ] += gCoord[ J_AXIS ];
        rightHandSide[ K_AXIS ] += gCoord[ K_AXIS ];

        /* Solve for xi increment */
        TensorArray_SolveSystem( jacobiMatrix, xiIncrement, rightHandSide, MT_VOLUME );

        /* Update xi */
        lCoord[ I_AXIS ] += xiIncrement[ I_AXIS ];
	lCoord[ J_AXIS ] += xiIncrement[ J_AXIS ];
        lCoord[ K_AXIS ] += xiIncrement[ K_AXIS ];

        /* Check for convergence */
        maxResidual = fabs( xiIncrement[ I_AXIS ] );
        if( maxResidual < fabs( xiIncrement[ J_AXIS ] ) )
            maxResidual = fabs( xiIncrement[ J_AXIS ] );
        if( maxResidual < fabs( xiIncrement[ K_AXIS ] ) )
            maxResidual = fabs( xiIncrement[ K_AXIS ] );

        if( maxResidual < tolerance )
            return;
    }
    /* if we are here, it means the iterative method didn't converge.
       Thus we set the local coord's to be invalid, i.e. greater than 1.0 */
    lCoord[0] = 1.1;
}

void SLIntegrator_Spherical_InitPatches( void* slIntegrator ) {
    SLIntegrator_Spherical* 	self 		= (SLIntegrator_Spherical*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    Grid*			vertGrid	= ((CartesianGenerator*)feMesh->generator)->vertGrid;
    unsigned			el_i, IJK[3], IJK_test[3], lNode_i, gNode_i, index;
    int				*inc, nInc;

    for( el_i = 0; el_i < Mesh_GetDomainSize( feMesh, MT_VOLUME ); el_i++ ) {
        FeMesh_GetElementNodes( feMesh, el_i, self->velocityField->inc );
        nInc = IArray_GetSize( self->velocityField->inc );
        inc  = IArray_GetPtr( self->velocityField->inc );

        gNode_i = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
        Grid_Lift( vertGrid, gNode_i, IJK );

        if( !SemiLagrangianIntegrator_HasRight3D( feMesh, self->inc, el_i, inc, nInc ) ) IJK[0]--;
        if( !SemiLagrangianIntegrator_HasTop3D(   feMesh, self->inc, el_i, inc, nInc ) ) IJK[1]--;
        if( !SemiLagrangianIntegrator_HasBack3D(  feMesh, self->inc, el_i, inc, nInc ) ) IJK[2]--;

        if( SemiLagrangianIntegrator_HasLeft3D(   feMesh, self->inc, el_i, inc, nInc ) ) IJK[0]--;
        if( SemiLagrangianIntegrator_HasBottom3D( feMesh, self->inc, el_i, inc, nInc ) ) IJK[1]--;
        if( SemiLagrangianIntegrator_HasFront3D(  feMesh, self->inc, el_i, inc, nInc ) ) IJK[2]--;
        FeMesh_NodeGlobalToDomain( feMesh, Grid_Project( vertGrid, IJK ), &lNode_i );

        index = 0;
        /* interpolate using Lagrange polynomial shape functions */
        for( IJK_test[2] = IJK[2]; IJK_test[2] < IJK[2] + 4; IJK_test[2]++ ) {
            for( IJK_test[1] = IJK[1]; IJK_test[1] < IJK[1] + 4; IJK_test[1]++ ) {
                for( IJK_test[0] = IJK[0]; IJK_test[0] < IJK[0] + 4; IJK_test[0]++ ) {
                    gNode_i = Grid_Project( vertGrid, IJK_test );
                    if( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_i, &lNode_i ) ) abort();
                    else
                        self->elPatch[el_i][index++] = lNode_i;
                }
            }
        }
    }
}

void SLIntegrator_Spherical_InitPatches_Quad( void* slIntegrator ) {
    SLIntegrator_Spherical* 	self 		= (SLIntegrator_Spherical*)slIntegrator;
    FeMesh*			feMesh		= self->velocityField->feMesh;
    Grid*                       grid            = ((CartesianGenerator*)feMesh->generator)->vertGrid;
    IArray*                     iArray1         = IArray_New();
    IArray*                     iArray2         = IArray_New();
    unsigned                    el_i, lNode, gNode;
    unsigned                    IJK[3], IJK_test[3], index;
    int                         nInc, *inc;
    int                         gOrig[8], orig_i;

    for( el_i = 0; el_i < Mesh_GetDomainSize( feMesh, MT_FACE ); el_i++ ) {
        FeMesh_GetElementNodes( feMesh, el_i, iArray1 );
        nInc = IArray_GetSize( iArray1 );
        inc  = IArray_GetPtr( iArray1 );

        gNode = Mesh_DomainToGlobal( feMesh, MT_VERTEX, inc[0] );
        /* left-bottom-front */
        Grid_Lift( grid, gNode, IJK );
        if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[0] = Grid_Project( grid, IJK );
        /* right-bottom-front */
        Grid_Lift( grid, gNode, IJK );
        if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[1] = Grid_Project( grid, IJK );
        /* left-top-front */
        Grid_Lift( grid, gNode, IJK );
        if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[2] = Grid_Project( grid, IJK );
        /* right-top-front */
        Grid_Lift( grid, gNode, IJK );
        if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if(  SemiLagrangianIntegrator_HasFront3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[3] = Grid_Project( grid, IJK );
        /* left-bottom-back */
        Grid_Lift( grid, gNode, IJK );
        if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[4] = Grid_Project( grid, IJK );
        /* right-bottom-back */
        Grid_Lift( grid, gNode, IJK );
        if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if(  SemiLagrangianIntegrator_HasBottom3D( feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[5] = Grid_Project( grid, IJK );
        /* left-top-back */
        Grid_Lift( grid, gNode, IJK );
        if(  SemiLagrangianIntegrator_HasLeft3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[6] = Grid_Project( grid, IJK );
        /* right-top-back */
        Grid_Lift( grid, gNode, IJK );
        if( !SemiLagrangianIntegrator_HasRight3D(  feMesh, iArray2, el_i, inc, nInc ) ) IJK[0]--;
        if( !SemiLagrangianIntegrator_HasTop3D(    feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        if( !SemiLagrangianIntegrator_HasBack3D(   feMesh, iArray2, el_i, inc, nInc ) ) IJK[1]--;
        gOrig[7] = Grid_Project( grid, IJK );

        for( orig_i = 0; orig_i < 8; orig_i++ ) {
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
    Stg_Class_Delete( iArray1 );
    Stg_Class_Delete( iArray2 );
}

double SLIntegrator_Spherical_CalcAdvDiffDt( void* slIntegrator, FiniteElementContext* context ) {
    SLIntegrator_Spherical*	self 		= (SLIntegrator_Spherical*)slIntegrator;
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
