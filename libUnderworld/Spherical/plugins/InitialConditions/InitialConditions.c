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
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Spherical/Spherical.h>

#include "InitialConditions.h"

const Type Spherical_InitialConditions_Type = "Spherical_InitialConditions";
Spherical_InitialConditions* Spherical_InitialConditions_selfPointer = NULL;

typedef double (*funcPtr)( double* coord, double xo, double yo, double ax, double ay );

double SphericalICs_SolWave( double* coord, double xo, double yo, double ax, double ay ) {
    double	rs[3];

    Spherical_XYZ2regionalSphere( coord, rs );

    double	sx 	= 1.0/cosh( ax*(rs[1] - xo) );
    double	sy 	= 1.0/cosh( ay*(rs[2] - yo) );

    return sx*sx*sy*sy;
}

void InitialConditions_SolWave( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*          	feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double			xo   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftEta" , 0.5 );
    double			yo		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_shiftZeta", 0.5 );
    double			ax   		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleEta" , 5.0 );
    double			ay		= Dictionary_GetDouble_WithDefault( context->dictionary, "solWave_scaleZeta", 7.0 );

    *result = SphericalICs_SolWave( coord, xo, yo, ax, ay );
}

void InitialConditions_ShearCellEta( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*          	feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double			rs[3];

    Spherical_XYZ2regionalSphere( coord, rs );

    *result = M_PI*sin( 0.5*M_PI*rs[1] )*cos( 0.5*M_PI*rs[2] );
}

void InitialConditions_ShearCellZeta( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*          	feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double			rs[3];

    Spherical_XYZ2regionalSphere( coord, rs );

    *result = -M_PI*cos( 0.5*M_PI*rs[1] )*sin( 0.5*M_PI*rs[2] );
}

void InitialConditions_ParametricSphereX( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*                 feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    double                      scale           = Dictionary_GetDouble_WithDefault( context->dictionary, "parametricSphere_scale", 1.0 );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double                      radius          = sqrt( coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2] );
    double                      phi;

    phi = atan2( coord[1], coord[0] );

    *result = -radius*sin( phi )/scale;
    //if( coord[0] < 0.0 )
    //    *result *= -1.0;
}

void InitialConditions_ParametricSphereY( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*                 feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    double                      scale           = Dictionary_GetDouble_WithDefault( context->dictionary, "parametricSphere_scale", 1.0 );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double                      radius          = sqrt( coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2] );
    double                      phi;

    phi = atan2( coord[1], coord[0] );

    *result = radius*cos( phi )/scale;
    //if( coord[0] < 0.0 )
    //    *result *= -1.0;
}

void InitialConditions_ParametricSphereZ( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _data, void* _result ) {
    FiniteElementContext*       context         = (FiniteElementContext*)_context;
    FeVariable*                 feVariable      = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "VelocityField" );
    double                      scale           = Dictionary_GetDouble_WithDefault( context->dictionary, "parametricSphere_scale", 1.0 );
    FeMesh*                     mesh            = feVariable->feMesh;
    double*                     coord           = Mesh_GetVertex( mesh, node_lI );
    double*                     result          = (double*)_result;
    double                      radius          = sqrt( coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2] );
    double                      phi;

    phi = atan2( coord[1], coord[0] );

    *result = 0.0;
}

void Spherical_InitialConditions_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
    Spherical_InitialConditions* 	self 		= (Spherical_InitialConditions*)_self;
    AbstractContext*		context		= Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, NULL  );
    ConditionFunction*		condFunc;

    condFunc = ConditionFunction_New( InitialConditions_SolWave, (Name)"SolWave_RS", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( InitialConditions_ShearCellEta, (Name)"ShearCellEta", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( InitialConditions_ShearCellZeta, (Name)"ShearCellZeta", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( InitialConditions_ParametricSphereX, (Name)"ParametricSphereX", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( InitialConditions_ParametricSphereY, (Name)"ParametricSphereY", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
    condFunc = ConditionFunction_New( InitialConditions_ParametricSphereZ, (Name)"ParametricSphereZ", NULL  );
    ConditionFunction_Register_Add( condFunc_Register, condFunc );
}

void Spherical_InitialConditions_Build( void* _self, void* data ) {
    Spherical_InitialConditions* self = (Spherical_InitialConditions*)_self;
}

void Spherical_InitialConditions_Initialise( void* _self, void* data ) {
    Spherical_InitialConditions* 	self 	= (Spherical_InitialConditions*)_self;
}

void Spherical_InitialConditions_Execute( void* _self, void* data ) {
    Spherical_InitialConditions* self = (Spherical_InitialConditions*)_self;
}

void* Spherical_InitialConditions_New( Name name ) {
    /* Variables set in this function */
    SizeT                                              _sizeOfSelf = sizeof(Spherical_InitialConditions);
    Type                                                      type = Spherical_InitialConditions_Type;
    Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
    Stg_Class_PrintFunction*                                _print = _Codelet_Print;
    Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
    Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Spherical_InitialConditions_New;
    Stg_Component_ConstructFunction*                    _construct = Spherical_InitialConditions_AssignFromXML;
    Stg_Component_BuildFunction*                            _build = Spherical_InitialConditions_Build;
    Stg_Component_InitialiseFunction*                  _initialise = Spherical_InitialConditions_Initialise;
    Stg_Component_ExecuteFunction*                        _execute = Spherical_InitialConditions_Execute;
    Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

    /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
    AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

    return _Codelet_New(  CODELET_PASSARGS  );
}

Index Spherical_InitialConditions_Register( PluginsManager* mgr ) {
    return PluginsManager_Submit( mgr, Spherical_InitialConditions_Type, (Name)"0", Spherical_InitialConditions_New );
}

