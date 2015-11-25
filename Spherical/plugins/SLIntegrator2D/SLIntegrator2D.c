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

#include "SLIntegrator2D.h"

const Type Spherical_SLIntegrator2D_Type = "Spherical_SLIntegrator2D";
Spherical_SLIntegrator2D* Spherical_SLIntegrator2D_selfPointer = NULL;

typedef double (*funcPtr)( double* coord, double alpha, double time );

double AngleVaryingTemp( double* coord, double alpha, double time ) {
    double	theta	= atan2( coord[1], coord[0] );

    return cos( alpha*theta );
}

double SLIntegrator2D_Dt( void* _context ) {
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

double SLIntegrator2D_EvaluateError( FiniteElementContext* context, FeVariable* phiField, Swarm* gaussSwarm, funcPtr func ) {
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
    double			alpha		= Dictionary_GetDouble_WithDefault( context->dictionary, "angleVaryingTemperature_alpha", 7.0 );
    double			time		= context->timeStep*SLIntegrator2D_Dt( context );
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

            initialValue = func( gCoord, alpha, time );
            
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

void SLIntegrator2D_UpdatePositions( void* data, FiniteElementContext* context ) {
    Swarm*			gaussSwarm      = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"gaussSwarm" );
    FeVariable*          	phiField	= (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureField" );
    double			error;
    FILE* 			fp 		= fopen( "error.txt", "w" );

    if( context->timeStep == context->maxTimeSteps ) {
        error = SLIntegrator2D_EvaluateError( context, phiField, gaussSwarm, AngleVaryingTemp );
        printf( "\nstep: %u\terror: %12.10f\n\n", context->timeStep, error );

        fp = fopen( "error.txt", "w" );
        fprintf( fp, "%12.10e", error );
        fclose( fp );
    }
}

void Spherical_SLIntegrator2D_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    Spherical_SLIntegrator2D* 	self 		= (Spherical_SLIntegrator2D*)_self;

    self->phiField = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureField" );
    self->context  = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, NULL  );

    ContextEP_ReplaceAll( self->context, AbstractContext_EP_Dt, SLIntegrator2D_Dt );
    ContextEP_Append( self->context, AbstractContext_EP_UpdateClass, SLIntegrator2D_UpdatePositions );
}

void Spherical_SLIntegrator2D_Build( void* _self, void* data ) {}

void Spherical_SLIntegrator2D_Initialise( void* _self, void* data ) {}

void Spherical_SLIntegrator2D_Execute( void* _self, void* data ) {}

void* Spherical_SLIntegrator2D_New( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Spherical_SLIntegrator2D);
    Type                                                      type = Spherical_SLIntegrator2D_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Spherical_SLIntegrator2D_New;
    Stg_Component_ConstructFunction*                    _construct = Spherical_SLIntegrator2D_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = Spherical_SLIntegrator2D_Build;
    Stg_Component_InitialiseFunction*                  _initialise = Spherical_SLIntegrator2D_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = Spherical_SLIntegrator2D_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

    return _Codelet_New(  CODELET_PASSARGS  );
}

Index Spherical_SLIntegrator2D_Register( PluginsManager* mgr ) {
    return PluginsManager_Submit( mgr, Spherical_SLIntegrator2D_Type, (Name)"0", Spherical_SLIntegrator2D_New );
}

