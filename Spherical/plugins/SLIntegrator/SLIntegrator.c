/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC)
** Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <Spherical/Spherical.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Spherical/Spherical.h>

#include "SLIntegrator.h"

const Type Spherical_SLIntegrator_Type = "Spherical_SLIntegrator";
Spherical_SLIntegrator* Spherical_SLIntegrator_selfPointer = NULL;

typedef double (*funcPtr)( double* coord, double xo, double yo, double ax, double ay );

double SolWave_Domain( double* coord, double xo, double yo, double ax, double ay ) {
    //double      rs[3];
    //Spherical_XYZ2regionalSphere( coord, rs );

    double      sx      = 1.0/cosh( ax*(coord[0] - xo) );
    double      sy      = 1.0/cosh( ay*(coord[1] - yo) );

    return sx*sx*sy*sy;
}

double SLIntegrator_Dt( void* _context ) {
    FiniteElementContext*       context         = (FiniteElementContext*) _context;
    double                      dt;
    FeVariable*          	velocityField   = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
    double			timeStepScale   = Dictionary_GetDouble_WithDefault( context->dictionary, "timeStepScale", 1.0 );
    double			staticTimeStep  = Dictionary_GetDouble_WithDefault( context->dictionary, "staticTimeStep", -1.0 );
    unsigned			circleStep      = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "circleStep", 0 );
    unsigned 			nSteps 		= Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "maxTimeSteps", 1 );
    double                      velMax;
    double                      delta[3], minDelta;
    double			dtMax;

    if( staticTimeStep > 0.0 ) {
        dt = staticTimeStep;
    }
    else {
        velMax = FieldVariable_GetMaxGlobalFieldMagnitude( velocityField );
        FeVariable_GetMinimumSeparation( velocityField, &minDelta, delta );

        dt = 0.5 * minDelta / velMax;
        dt *= timeStepScale;
    }

    if( circleStep ) {
        dt    = 2.0*M_PI*0.25/nSteps;
        dtMax = delta[0]/velMax;
        if( dt > dtMax ) {
            printf( "ERROR - dt: %12.10e > dt_max: %12.10e\n", dt, dtMax );
            abort();
        }
    }

    printf( "time step: %12.10f\n", dt );

    return dt;
}

double SLIntegrator_EvaluateError( FiniteElementContext* context, FeVariable* phiField, Swarm* gaussSwarm, funcPtr func ) {
    FeMesh*              	feMesh          = phiField->feMesh;
    GaussParticleLayout* 	particleLayout  = (GaussParticleLayout*)gaussSwarm->particleLayout;
    Index                	lElement_I, lCell_I;
    unsigned             	nDims           = Mesh_GetDimSize( feMesh );
    unsigned             	numMeshElements = Mesh_GetLocalSize( feMesh, nDims );
    double               	lError          = 0.0;
    double               	lAnalytic       = 0.0;
    double               	gError, gAnalytic, gErrorNorm;
    IntegrationPoint*    	gaussPoint;
    unsigned             	gaussPoint_I, numGaussPoints;
    double               	initialValue, finalValue;
    double               	elErrorSq, elAnalyticSq;
    ElementType*         	elementType;
    double               	detJac;
    double			gCoord[3];
    double			xo   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftX", 0.5 );
    double			yo		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftY", 0.5 );
    double			ax   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleX", 5.0 );
    double			ay		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleY", 7.0 );
    void*     			slIntegrator    = LiveComponentRegister_Get( context->CF->LCRegister, (Name)"integrator" );

    for( lElement_I = 0; lElement_I < numMeshElements; lElement_I++ ) {
        lCell_I = CellLayout_MapElementIdToCellId( gaussSwarm->cellLayout, lElement_I );
        numGaussPoints = _GaussParticleLayout_InitialCount( particleLayout, NULL, lCell_I );

        elementType = FeMesh_GetElementType( feMesh, lElement_I );

        elErrorSq    = 0.0;
        elAnalyticSq = 0.0;

        for( gaussPoint_I = 0; gaussPoint_I < numGaussPoints; gaussPoint_I++ ) {
            gaussPoint = (IntegrationPoint*) Swarm_ParticleInCellAt( gaussSwarm, lCell_I, gaussPoint_I );
            FeMesh_CoordLocalToGlobal( feMesh, lElement_I, gaussPoint->xi, gCoord );

            initialValue = func( gCoord, xo, yo, ax, ay );
            
            if( nDims == 3 )
                SLIntegrator_Spherical_CubicInterpolator( slIntegrator, phiField, gCoord, &finalValue );
            else
                SLIntegrator_Polar_CubicInterpolator( slIntegrator, phiField, gCoord, &finalValue );
 
            detJac = ElementType_JacobianDeterminant( elementType, feMesh, lElement_I, gaussPoint->xi, nDims );

            elErrorSq    += ( finalValue - initialValue ) * ( finalValue - initialValue ) * gaussPoint->weight * detJac;
            elAnalyticSq += ( initialValue * initialValue ) * gaussPoint->weight * detJac;
        }

        lError    += sqrt( elErrorSq );
        lAnalytic += sqrt( elAnalyticSq );
    }

    MPI_Allreduce( &lError,    &gError,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( &lAnalytic, &gAnalytic, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    gErrorNorm = gError / gAnalytic;

    return gErrorNorm;
}

void SLIntegrator_UpdatePositions( void* data, FiniteElementContext* context ) {
    Index			reverseTimeStep = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "reverseTimeStep", 100 );
    FeVariable*          	velocityField   = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
    Swarm*			gaussSwarm      = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"gaussSwarm" );
    FeMesh*                     mesh            = velocityField->feMesh;
    unsigned             	node_I;
    Index                       dim_I;
    unsigned             	nDims           = Mesh_GetDimSize( mesh );
    double                      velocity[3];
    double                      phi;
    FeVariable*          	phiField	= (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureField" );
    FeVariable*          	phiInitField	= NULL;
    double			error;
    double			xo   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftX", 0.5 );
    double			yo		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftY", 0.5 );
    double			ax   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleX", 5.0 );
    double			ay		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleY", 7.0 );

    _FeVariable_SyncShadowValues( velocityField );

    /* reverse the numerically advected particles (& the semi lagrangian field also) */
    if( context->timeStep == reverseTimeStep ) {
        for( node_I = 0; node_I < Mesh_GetLocalSize( mesh, MT_VERTEX ); node_I++ ) {
            _FeVariable_GetValueAtNode( velocityField, node_I, velocity );
 
           for( dim_I = 0; dim_I < nDims; dim_I++ ) {
               velocity[dim_I] *= -1;
           }

           FeVariable_SetValueAtNode( velocityField, node_I, velocity );
        }
    }

    if( context->timeStep == context->maxTimeSteps ) {
        phiInitField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureInitField" );
        error        = SLIntegrator_EvaluateError( context, phiField, gaussSwarm, SolWave_Domain );
        printf( "\nstep: %u\terror: %12.10f\n\n", context->timeStep, error );

        FILE* fp = fopen( "error.txt", "w" );
        fprintf( fp, "%12.10e", error );
        fclose( fp );

        double phi_a;
        for( node_I = 0; node_I < Mesh_GetLocalSize( mesh, MT_VERTEX ); node_I++ ) {
            FeVariable_GetValueAtNode( phiField, node_I, &phi );
            phi_a = SolWave_Domain( Mesh_GetVertex( mesh, node_I ), xo, yo, ax, ay );
            phi_a -= phi;
            FeVariable_SetValueAtNode( phiInitField, node_I, &phi_a );
        }
    }
}

void Spherical_SLIntegrator_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    Spherical_SLIntegrator* 	self 		= (Spherical_SLIntegrator*)_self;
    AbstractContext*		context		= Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, NULL  );

    ContextEP_ReplaceAll( context, AbstractContext_EP_Dt, SLIntegrator_Dt );
    ContextEP_Append( context, AbstractContext_EP_UpdateClass, SLIntegrator_UpdatePositions );

    self->phiField     = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureField" );
    self->phiInitField = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureInitField" );
}

void Spherical_SLIntegrator_Build( void* _self, void* data ) {
    Spherical_SLIntegrator* self = (Spherical_SLIntegrator*)_self;
}

void Spherical_SLIntegrator_Initialise( void* _self, void* data ) {
    Spherical_SLIntegrator* 	self 	= (Spherical_SLIntegrator*)_self;
    double 			phi[3];
    int				node_i;

    for( node_i = 0; node_i < Mesh_GetLocalSize( self->phiField->feMesh, MT_VERTEX ); node_i++ ) {
        FeVariable_GetValueAtNode( self->phiField, node_i, phi );
        FeVariable_SetValueAtNode( self->phiInitField, node_i, phi );
    }
}

void Spherical_SLIntegrator_Execute( void* _self, void* data ) {
    Spherical_SLIntegrator* self = (Spherical_SLIntegrator*)_self;

}

void* Spherical_SLIntegrator_New( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Spherical_SLIntegrator);
    Type                                                      type = Spherical_SLIntegrator_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Spherical_SLIntegrator_New;
    Stg_Component_ConstructFunction*                    _construct = Spherical_SLIntegrator_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = Spherical_SLIntegrator_Build;
    Stg_Component_InitialiseFunction*                  _initialise = Spherical_SLIntegrator_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = Spherical_SLIntegrator_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

    return _Codelet_New(  CODELET_PASSARGS  );
}

Index Spherical_SLIntegrator_Register( PluginsManager* mgr ) {
    return PluginsManager_Submit( mgr, Spherical_SLIntegrator_Type, (Name)"0", Spherical_SLIntegrator_New );
}

