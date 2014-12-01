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
#include <StgFEM/StgFEM.h>

#include "SemiLagrangianADE_CalcErrors.h"

const Type StgFEM_SemiLagrangianADE_CalcErrors_Type = "StgFEM_SemiLagrangianADE_CalcErrors";
StgFEM_SemiLagrangianADE_CalcErrors* StgFEM_SemiLagrangianADE_CalcErrors_selfPointer = NULL;

typedef double (*funcPtr)( FiniteElementContext* context, double* coord );

double SemiLagrangianADE_CalcErrors_Dt( void* _context ) {
    FiniteElementContext*       context         = (FiniteElementContext*) _context;
    FeVariable*                 velocityField   = (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
    double                      courantNumber   = Dictionary_GetDouble_WithDefault( context->dictionary, "courantNumber", 1.0 );
    double                      velMax          = Dictionary_GetDouble_WithDefault( context->dictionary, "velMax", 1.0 );
    double			manualDt	= Dictionary_GetDouble_WithDefault( context->dictionary, "manualTimeStep", 0.0 );
    double			velMaxObs;
    double                      delta[3], minDelta;
    double                      dt, dtMax;

    if( manualDt > 1.0e-6 ) {
        return manualDt;
    }

    velMaxObs = FieldVariable_GetMaxGlobalFieldMagnitude( velocityField );
    FeVariable_GetMinimumSeparation( velocityField, &minDelta, delta );

    dt    = courantNumber*minDelta/velMax;
    dtMax = minDelta/velMaxObs;

    if( dt > dtMax ) {
        printf( "\nERROR: time step, %12.10e > maximum allowed, %12.10e\n\n", dt, dtMax );
    }

    return dt;
}

void SemiLagrangianADE_CalcErrors_ShearCellX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*          	feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;

    *result = M_PI * sin( M_PI * coord[0] ) * cos( M_PI * coord[1] );
}

void SemiLagrangianADE_CalcErrors_ShearCellY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*          	feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;

    *result = -M_PI * cos( M_PI * coord[0] ) * sin( M_PI * coord[1] );
}

double FourierMode( FiniteElementContext* context, double* coord ) {
        Dictionary    *dictionary = context->dictionary;
        FeVariable    *feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
	double        time        = context->timeStep*SemiLagrangianADE_CalcErrors_Dt( context );
	//double        time        = context->currentTime;
        double        kx, ky, vel[2], visc;
        double	      gCoord[3];
        double	      delta[3];
        unsigned      nNodes[3];

        kx   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"k_x",  1.0 );
        ky   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"k_y",  1.0 );
        visc = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"visc", 1.0 );

	kx *= 2.0*M_PI;
	ky *= 2.0*M_PI;

        SemiLagrangianIntegrator_GetDeltaConst( feVariable, delta, nNodes );
        SemiLagrangianIntegrator_BicubicInterpolator( feVariable, gCoord, delta, nNodes, vel );

        return exp( -visc*(kx*kx + ky*ky)*time )*cos( kx*coord[0] + ky*coord[1] - (vel[0]*kx + vel[1]*ky)*time );
}

void SemiLagrangianADE_CalcErrors_FourierMode(Node_LocalIndex node_lI,Variable_Index var_I,void *_context,void* _data, void* _result) {
        FiniteElementContext	*context    = (FiniteElementContext*)_context;
        Dictionary    		*dictionary = context->dictionary;
        FeVariable    		*feVariable = NULL;
        FeMesh        		*mesh       = NULL;
        double        		*result     = (double*) _result;
        double        		*coord;
        double        		kx, ky, vel[2], visc, time = 0.0;

        feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
        mesh       = feVariable->feMesh;

        kx   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"k_x",  1.0 );
        ky   = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"k_y",  1.0 );
        visc = Dictionary_GetDouble_WithDefault( dictionary, (Dictionary_Entry_Key)"visc", 1.0 );

	kx *= 2.0*M_PI;
	ky *= 2.0*M_PI;

        FeVariable_GetValueAtNode( feVariable, node_lI, vel );

        /* Find coordinate of node */
        coord = Mesh_GetVertex( mesh, node_lI );

        *result = exp( -visc*(kx*kx + ky*ky)*time )*cos( kx*coord[0] + ky*coord[1] - (vel[0]*kx + vel[1]*ky)*time );
}

double SemiLagrangianADE_CalcErrors_EvaluateError( FiniteElementContext* context, FeVariable* phiField, Swarm* gaussSwarm, funcPtr func ) {
    FeMesh*              	feMesh          = phiField->feMesh;
    GaussParticleLayout* 	particleLayout  = (GaussParticleLayout*)gaussSwarm->particleLayout;
    Index                	lElement_I, lCell_I;
    unsigned			nDims		= Mesh_GetDimSize( feMesh );
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
    double			delta[3];
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

            initialValue = func( context, gCoord );
            
            SemiLagrangianIntegrator_GetDeltaConst( phiField, delta, nNodes );
            SemiLagrangianIntegrator_BicubicInterpolator( phiField, gCoord, delta, nNodes, &finalValue );
 
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

void SemiLagrangianADE_CalcErrors_Evaluate( void* data, FiniteElementContext* context ) {
    FeVariable*          	velocityField   	= (FeVariable*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"VelocityField" );
    Swarm*			gaussSwarm      	= (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"gaussSwarm" );
    FeMesh*                     mesh            	= velocityField->feMesh;
    unsigned             	node_I;
    FeVariable*          	temperatureField	= NULL;
    FeVariable*          	tempFinalField		= NULL;
    double			error;
    double			temp_a, temp_n, diff;

    _FeVariable_SyncShadowValues( velocityField );

    temperatureField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TemperatureField" );
    tempFinalField   = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)"TempFinalField" );

    printf( "\ndt: %12.10e\n\n", SemiLagrangianADE_CalcErrors_Dt( context ) );

    for( node_I = 0; node_I < Mesh_GetLocalSize( mesh, MT_VERTEX ); node_I++ ) {
        temp_a = FourierMode( context, Mesh_GetVertex( temperatureField->feMesh, node_I ) );
        FeVariable_GetValueAtNode( temperatureField, node_I, &temp_n );
        diff = temp_a - temp_n;
        FeVariable_SetValueAtNode( tempFinalField, node_I, &temp_a );
    }

    if( context->timeStep == context->maxTimeSteps ) {
        error          = SemiLagrangianADE_CalcErrors_EvaluateError( context, temperatureField, gaussSwarm, FourierMode );
        printf( "\nstep: %u\terror: %12.10f\n\n", context->timeStep, error );

        FILE* fp = fopen( "error.txt", "w" );
        fprintf( fp, "%12.10e", error );
        fclose( fp );
    }
}

void StgFEM_SemiLagrangianADE_CalcErrors_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    StgFEM_SemiLagrangianADE_CalcErrors* 		self 		= (StgFEM_SemiLagrangianADE_CalcErrors*)_self;
    AbstractContext*		context		= Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, NULL  );
    ConditionFunction*		condFunc;

    condFunc = ConditionFunction_New( SemiLagrangianADE_CalcErrors_FourierMode, (Name)"FourierMode", NULL );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( SemiLagrangianADE_CalcErrors_ShearCellX, (Name)"ShearCellX", NULL );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( SemiLagrangianADE_CalcErrors_ShearCellY, (Name)"ShearCellY", NULL );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );

    ContextEP_ReplaceAll( context, AbstractContext_EP_Dt, SemiLagrangianADE_CalcErrors_Dt );
    ContextEP_Append( context, AbstractContext_EP_UpdateClass, SemiLagrangianADE_CalcErrors_Evaluate );
}

void StgFEM_SemiLagrangianADE_CalcErrors_Build( void* _self, void* data ) {
    StgFEM_SemiLagrangianADE_CalcErrors* self = (StgFEM_SemiLagrangianADE_CalcErrors*)_self;
}

void StgFEM_SemiLagrangianADE_CalcErrors_Initialise( void* _self, void* data ) {
    StgFEM_SemiLagrangianADE_CalcErrors* 	self 	= (StgFEM_SemiLagrangianADE_CalcErrors*)_self;
}

void StgFEM_SemiLagrangianADE_CalcErrors_Execute( void* _self, void* data ) {
    StgFEM_SemiLagrangianADE_CalcErrors* self = (StgFEM_SemiLagrangianADE_CalcErrors*)_self;

}

void* StgFEM_SemiLagrangianADE_CalcErrors_New( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(StgFEM_SemiLagrangianADE_CalcErrors);
    Type                                                      type = StgFEM_SemiLagrangianADE_CalcErrors_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = StgFEM_SemiLagrangianADE_CalcErrors_New;
    Stg_Component_ConstructFunction*                    _construct = StgFEM_SemiLagrangianADE_CalcErrors_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = StgFEM_SemiLagrangianADE_CalcErrors_Build;
    Stg_Component_InitialiseFunction*                  _initialise = StgFEM_SemiLagrangianADE_CalcErrors_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = StgFEM_SemiLagrangianADE_CalcErrors_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

    return _Codelet_New(  CODELET_PASSARGS  );
}

Index StgFEM_SemiLagrangianADE_CalcErrors_Register( PluginsManager* mgr ) {
    return PluginsManager_Submit( mgr, StgFEM_SemiLagrangianADE_CalcErrors_Type, (Name)"0", StgFEM_SemiLagrangianADE_CalcErrors_New );
}

