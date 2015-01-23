/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the SemiLagrangianIntegratorSuite
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "SemiLagrangianIntegratorSuite.h"

#define CURR_MODULE_NAME "SemiLagrangianIntegratorSuite"

typedef struct {
} SemiLagrangianIntegratorSuiteData;

typedef double (*funcPtr)( double* coord );

double SolWave( double* coord ) {
    double      sx      = 1.0/cosh( 5.0*(coord[0] - 0.5) );
    double      sy      = 1.0/cosh( 7.0*(coord[1] - 0.5) );

    return sx*sx*sy*sy;
}

double SemiLagrangianIntegratorSuite_Dt( void* _context ) {
    FiniteElementContext*	context		= (FiniteElementContext*) _context;
    FeVariable*			velocityField	= (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField"  );
    Index                   	staticTimeStep  = Dictionary_GetDouble_WithDefault( context->dictionary, "staticTimeStep", 0.0 );
    double			velMax;
    double			delta[3], minDelta;

    if( staticTimeStep > 1.0e-6 )
        return staticTimeStep;

    velMax = FieldVariable_GetMaxGlobalFieldMagnitude( velocityField );
    //FeVariable_GetMinimumSeparation( velocityField, &minDelta, delta );
    minDelta = 1.0/512;

    return 0.5 * minDelta / velMax;
}

void SemiLagrangianIntegratorSuite_SolWave( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*                 feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;

    *result = SolWave( coord );
}

void SemiLagrangianIntegratorSuite_ShearCellX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*	context		= (FiniteElementContext*)_context;
    FeVariable*			feVariable	= (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*			mesh		= feVariable->feMesh;
    double*			coord		= Mesh_GetVertex( mesh, node_lI );
    double*			result		= (double*)_result;

    *result = M_PI * sin( M_PI * coord[0] ) * cos( M_PI * coord[1] );
}

void SemiLagrangianIntegratorSuite_ShearCellY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*	context		= (FiniteElementContext*)_context;
    FeVariable*			feVariable	= (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*			mesh		= feVariable->feMesh;
    double*			coord		= Mesh_GetVertex( mesh, node_lI );
    double*			result		= (double*)_result;

    *result = -M_PI * cos( M_PI * coord[0] ) * sin( M_PI * coord[1] );
}

void SemiLagrangianIntegratorSuite_ParametricCircleX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*	context		= (FiniteElementContext*)_context;
    FeVariable*			feVariable	= (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*			mesh		= feVariable->feMesh;
    double*			coord		= Mesh_GetVertex( mesh, node_lI );
    double*			result		= (double*)_result;
    double			radius		= sqrt( coord[0]*coord[0] + coord[1]*coord[1] );
    double			theta		= atan( coord[1]/coord[0] );

    *result = -radius*sin( theta );
}

void SemiLagrangianIntegratorSuite_ParametricCircleY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*	context		= (FiniteElementContext*)_context;
    FeVariable*			feVariable	= (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*			mesh		= feVariable->feMesh;
    double*			coord		= Mesh_GetVertex( mesh, node_lI );
    double*			result		= (double*)_result;
    double			radius		= sqrt( coord[0]*coord[0] + coord[1]*coord[1] );
    double			theta		= atan( coord[1]/coord[0] );

    *result = radius*cos( theta );
}

void SemiLagrangianIntegratorSuite_UpdatePositions( void* data, FiniteElementContext* context ) {
    Index                   	reverseTimeStep = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "reverseTimeStep", 100 );
    FeVariable*			velocityField	= (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
    FeMesh*			mesh		= velocityField->feMesh;
    unsigned			node_I, dim_I;
    unsigned			nDims		= Mesh_GetDimSize( mesh );
    double			velocity[3];

    _FeVariable_SyncShadowValues( velocityField );

    /* reverse the numerically advected particles (& the semi lagrangian field also) */
    if( context->timeStep == reverseTimeStep ) {
        for( node_I = 0; node_I < Mesh_GetLocalSize( mesh, MT_VERTEX ); node_I++ ) {
            _FeVariable_GetValueAtNode( velocityField, node_I, velocity );

            for( dim_I = 0; dim_I < nDims; dim_I++ ) {
                velocity[dim_I] *= -1.0;
            }
            FeVariable_SetValueAtNode( velocityField, node_I, velocity );
        }
    }
    _FeVariable_SyncShadowValues( velocityField );
}

double SemiLagrangianIntegratorSuite_EvaluateError( SemiLagrangianIntegrator* slIntegrator, FeVariable* temperatureField, Swarm* gaussSwarm, funcPtr func ) {
    FeMesh*			feMesh		= temperatureField->feMesh;
    GaussParticleLayout*	particleLayout 	= (GaussParticleLayout*)gaussSwarm->particleLayout;
    Index			lElement_I, lCell_I;
    unsigned			nDims		= Mesh_GetDimSize( feMesh );
    unsigned			numMeshElements	= Mesh_GetLocalSize( feMesh, nDims );
    double			lError		= 0.0;
    double			lAnalytic 	= 0.0;
    double			gError, gAnalytic;
    IntegrationPoint*		gaussPoint;
    unsigned			gaussPoint_I, numGaussPoints;
    double			initialValue, finalValue;
    double			elErrorSq, elAnalyticSq;
    ElementType*		elementType;
    double			detJac;
    double			gCoord[3], delta[3];
    unsigned			nNodes[3];

    for( lElement_I = 0; lElement_I < numMeshElements; lElement_I++ ) {
        lCell_I = CellLayout_MapElementIdToCellId( gaussSwarm->cellLayout, lElement_I );
        numGaussPoints = _GaussParticleLayout_InitialCount( particleLayout, NULL, lCell_I );

        elementType = FeMesh_GetElementType( feMesh, lElement_I );

        elErrorSq    = 0.0;
        elAnalyticSq = 0.0;

        for( gaussPoint_I = 0; gaussPoint_I < numGaussPoints; gaussPoint_I++ ) {
            gaussPoint = (IntegrationPoint*) Swarm_ParticleInCellAt( gaussSwarm, lCell_I, gaussPoint_I );

            FeMesh_CoordLocalToGlobal( feMesh, lElement_I, gaussPoint->xi, gCoord );
            initialValue = func( gCoord );

            SemiLagrangianIntegrator_GetDeltaConst( temperatureField, delta, nNodes );
            SemiLagrangianIntegrator_CubicInterpolator( slIntegrator, temperatureField, gCoord, delta, nNodes, &finalValue );

            detJac = ElementType_JacobianDeterminant( elementType, feMesh, lElement_I, gaussPoint->xi, nDims );

            elErrorSq += ( finalValue - initialValue ) * ( finalValue - initialValue ) * gaussPoint->weight * detJac;
            elAnalyticSq += ( initialValue * initialValue ) * gaussPoint->weight * detJac;
        }

        lError    += sqrt( elErrorSq );
        lAnalytic += sqrt( elAnalyticSq );
    }

    MPI_Allreduce( &lError, &gError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( &lAnalytic, &gAnalytic, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    return gError/gAnalytic;
}

void SemiLagrangianIntegratorSuite_Setup( SemiLagrangianIntegratorSuiteData* data ) {}

void SemiLagrangianIntegratorSuite_Teardown( SemiLagrangianIntegratorSuiteData* data ) {}

void SemiLagrangianIntegratorSuite_Test( SemiLagrangianIntegratorSuiteData* data ) {
    Stg_ComponentFactory*	cf;
    ConditionFunction*      	condFunc;
    //char			xml_input[PCU_PATH_MAX];
    double			l2Error;
    FeVariable*			temperatureField;
    FeVariable*			temperatureInitField;
    Swarm*			gaussSwarm;
    double			temperature[3];
    unsigned			node_i;
    AbstractContext*		context;
    SemiLagrangianIntegrator*	slIntegrator;

    //pcu_filename_input( "testSemiLagrangianIntegrator.xml", xml_input );
    cf = stgMainInitFromXML( "StgFEM/Utils/input/testSemiLagrangianIntegrator.xml", MPI_COMM_WORLD, NULL );
    context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, NULL  );

    condFunc = ConditionFunction_New( SemiLagrangianIntegratorSuite_SolWave, (Name)"SolWave", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( SemiLagrangianIntegratorSuite_ShearCellX, (Name)"ShearCellX", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( SemiLagrangianIntegratorSuite_ShearCellY, (Name)"ShearCellY", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );

    /* manually set the timestep */
    ContextEP_ReplaceAll( context, AbstractContext_EP_Dt, SemiLagrangianIntegratorSuite_Dt );
    ContextEP_Append( context, AbstractContext_EP_UpdateClass, SemiLagrangianIntegratorSuite_UpdatePositions );

    stgMainBuildAndInitialise( cf );

    temperatureField     = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureField" );
    temperatureInitField = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureInitField" );
    gaussSwarm           = (Swarm*)LiveComponentRegister_Get( cf->LCRegister, (Name)"gaussSwarm" );
    slIntegrator         = (SemiLagrangianIntegrator*)LiveComponentRegister_Get( cf->LCRegister, (Name)"integrator" );
    for( node_i = 0; node_i < Mesh_GetLocalSize( temperatureField->feMesh, MT_VERTEX ); node_i++ ) {
        FeVariable_GetValueAtNode( temperatureField,     node_i, temperature );
        FeVariable_SetValueAtNode( temperatureInitField, node_i, temperature );
    }

    stgMainLoop( cf );

    l2Error = SemiLagrangianIntegratorSuite_EvaluateError( slIntegrator, temperatureField, gaussSwarm, SolWave );

    printf( "\ntime step: %12.10e\n\n", SemiLagrangianIntegratorSuite_Dt(context) );
    printf( "\nERROR (1): %12.10e\n\n", l2Error );

    pcu_check_true( l2Error < 5.0e-4 );

    stgMainDestroy( cf );
}

double Sech( double a, double x )    { return 1.0/cosh(a*(x-0.5)); }
double Sech2( double a, double x )   { return Sech(a,x)*Sech(a,x); }
double f( double a[2], double x[2] ) { return Sech2(a[0],x[0])*Sech2(a[1],x[1]); }

void SemiLagrangianIntegratorSuite_LagrangianInterpolation( SemiLagrangianIntegratorSuiteData* data ) {
    Stg_ComponentFactory*	cf;
    FeVariable*			temperatureField;
    Swarm*			gaussSwarm;
    ElementType*		elType;
    IntegrationPoint*		gaussPoint;
    Index			node_i, el_i, cell_i, pt_i;
    Index			nGaussPts;
    double			temperature, temperature_a;
    double			a[2]		= { 4.0, 6.0 };
    double			delta[2];
    unsigned			nNodes[2];
    double			gCoord[2];
    double			detJac, elErrorSq, elAnalyticSq;
    double			lError = 0.0, lAnalytic = 0.0, gError, gAnalytic, l2Error;
    SemiLagrangianIntegrator*	slIntegrator;

    cf = stgMainInitFromXML( "StgFEM/Utils/input/testSemiLagrangianIntegrator2.xml", MPI_COMM_WORLD, NULL );

    stgMainBuildAndInitialise( cf );

    slIntegrator     = (SemiLagrangianIntegrator*)LiveComponentRegister_Get( cf->LCRegister, (Name)"integrator" );
    temperatureField = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"TemperatureField" );
    gaussSwarm       = (Swarm*)LiveComponentRegister_Get( cf->LCRegister, (Name)"gaussSwarm" );

    SemiLagrangianIntegrator_GetDeltaConst( temperatureField, delta, nNodes );

    for( node_i = 0; node_i < Mesh_GetLocalSize( temperatureField->feMesh, MT_VERTEX ); node_i++ ) {
        temperature = f( a, Mesh_GetVertex( temperatureField->feMesh, node_i ) );
        FeVariable_SetValueAtNode( temperatureField, node_i, &temperature );
    }

    for( el_i = 0; el_i < Mesh_GetLocalSize( temperatureField->feMesh, MT_FACE ); el_i++ ) {
        cell_i    = CellLayout_MapElementIdToCellId( gaussSwarm->cellLayout, el_i );
        nGaussPts = _GaussParticleLayout_InitialCount( gaussSwarm->particleLayout, NULL, cell_i );
        elType    = FeMesh_GetElementType( temperatureField->feMesh, el_i );

        elErrorSq    = 0.0;
        elAnalyticSq = 0.0;

        for( pt_i = 0; pt_i < nGaussPts; pt_i++ ) {
            gaussPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( gaussSwarm, cell_i, pt_i );
            FeMesh_CoordLocalToGlobal( temperatureField->feMesh, el_i, gaussPoint->xi, gCoord );
            SemiLagrangianIntegrator_CubicInterpolator( slIntegrator, temperatureField, gCoord, delta, nNodes, &temperature );

            temperature_a = f( a, gCoord );
 
            detJac = ElementType_JacobianDeterminant( elType, temperatureField->feMesh, el_i, gaussPoint->xi, 2 );

            elErrorSq    += ( temperature - temperature_a ) * ( temperature - temperature_a ) * gaussPoint->weight * detJac;
            elAnalyticSq += ( temperature_a * temperature_a ) * gaussPoint->weight * detJac;
        }

        lError    += sqrt( elErrorSq );
        lAnalytic += sqrt( elAnalyticSq );
    }

    MPI_Allreduce( &lError,    &gError,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( &lAnalytic, &gAnalytic, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    l2Error = gError/gAnalytic;
    printf( "\nERROR (2): %12.10e\n\n", l2Error );

    pcu_check_true( l2Error < 1.0e-3 );
}

void SemiLagrangianIntegratorSuite_RK4Integration( SemiLagrangianIntegratorSuiteData* data ) {
    Stg_ComponentFactory*	cf;
    AbstractContext*		context;
    FeVariable*			velField;
    ConditionFunction*      	condFunc;
    unsigned                  	nSteps, step_i;
    double			dt;
    double			coord[2]	= { 0.75, 0.0 };
    double			cNew[2];
    int				pass;
    double			delta[2];
    unsigned			nnodes[2];
    SemiLagrangianIntegrator*	slIntegrator;

    cf = stgMainInitFromXML( "StgFEM/Utils/input/testSemiLagrangianIntegrator3.xml", MPI_COMM_WORLD, NULL );
    context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, NULL  );

    slIntegrator = (SemiLagrangianIntegrator*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"integrator" );

    condFunc = ConditionFunction_New( SemiLagrangianIntegratorSuite_ParametricCircleX, (Name)"ParametricCircleX", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( SemiLagrangianIntegratorSuite_ParametricCircleY, (Name)"ParametricCircleY", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );

    stgMainBuildAndInitialise( cf );

    velField = (FeVariable*)LiveComponentRegister_Get( cf->LCRegister, (Name)"VelocityField" );
    nSteps   = Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "maxTimeSteps", 1 );
    dt       = 0.5*M_PI/nSteps; //quarter circle

    SemiLagrangianIntegrator_GetDeltaConst( velField, delta, nnodes );

    for( step_i = 1; step_i <= nSteps; step_i++ ) {
        SemiLagrangianIntegrator_IntegrateRK4( slIntegrator, dt, delta, nnodes, coord, cNew );
        coord[0] = cNew[0];
        coord[1] = cNew[1];
    }

    printf( "\ndt:\t%12.10e\n", dt );
    printf( "\nposition: [%12.10e,%12.10e]\n\n", coord[0],coord[1] );
    
    pass = ( coord[0] - 0.0 < 1.0e-4 && coord[1] + 0.75 < 1.0e-4 ) ? 1 : 0;
    pcu_check_true( pass );
}

void SemiLagrangianIntegratorSuite( pcu_suite_t* suite ) {
    pcu_suite_setData( suite, SemiLagrangianIntegratorSuiteData );
    pcu_suite_setFixtures( suite, SemiLagrangianIntegratorSuite_Setup, SemiLagrangianIntegratorSuite_Teardown );
    pcu_suite_addTest( suite, SemiLagrangianIntegratorSuite_Test );
    pcu_suite_addTest( suite, SemiLagrangianIntegratorSuite_LagrangianInterpolation );
    pcu_suite_addTest( suite, SemiLagrangianIntegratorSuite_RK4Integration );
}

