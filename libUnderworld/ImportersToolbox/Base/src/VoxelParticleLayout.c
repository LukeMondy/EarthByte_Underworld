/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
** Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
** Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
** Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
** Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
** Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
** Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "VoxelDataHandler_Abstract.h"
#include "VoxelParticleLayout.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelParticleLayout_Type = "VoxelParticleLayout";

VoxelParticleLayout* _VoxelParticleLayout_New( VOXELPARTICLELAYOUT_DEFARGS  )
{
   VoxelParticleLayout* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelParticleLayout ) );
   self = (VoxelParticleLayout*)_GlobalParticleLayout_New(  GLOBALPARTICLELAYOUT_PASSARGS  );

   /** allocate small buffer memory.. for simplicity, choose size double to ensure enough memory */
   self->voxelData = Memory_Alloc_Bytes_Unnamed( sizeof(double), type );

   return self;
}

void _VoxelParticleLayout_Delete( void* particleLayout ) {
   VoxelParticleLayout* self = (VoxelParticleLayout*)particleLayout;

   Memory_Free(self->voxelData);

   /* Stg_Class_Delete parent class */
   _GlobalParticleLayout_Delete( self );
}

void _VoxelParticleLayout_Print( void* particleLayout, Stream* stream ) {
   VoxelParticleLayout* self = (VoxelParticleLayout*)particleLayout;

   Stream_Indent( stream );

   /* Parent class info */
   _GlobalParticleLayout_Print( self, stream );

   Stream_UnIndent( stream );
}

void* _VoxelParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
   VoxelParticleLayout*  self = (VoxelParticleLayout*)particleLayout;
   VoxelParticleLayout*  newVoxelParticleLayout;

   newVoxelParticleLayout = (VoxelParticleLayout*)_GlobalParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );

   return (void*)newVoxelParticleLayout;
}

void* _VoxelParticleLayout_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelParticleLayout);
   Type                                                                        type = VoxelParticleLayout_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelParticleLayout_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelParticleLayout_Print;
   Stg_Class_CopyFunction*                                                    _copy = _VoxelParticleLayout_Copy;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelParticleLayout_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelParticleLayout_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelParticleLayout_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelParticleLayout_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelParticleLayout_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelParticleLayout_Destroy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   ParticleLayout_SetInitialCountsFunction*                       _setInitialCounts = _VoxelParticleLayout_SetInitialCounts;
   ParticleLayout_InitialiseParticlesFunction*                 _initialiseParticles = _VoxelParticleLayout_InitialiseParticles;
   CoordSystem                                                          coordSystem = GlobalCoordSystem;
   Bool                                                 weightsInitialisedAtStartup = False;
   GlobalParticleLayout_InitialiseParticleFunction*             _initialiseParticle = _VoxelParticleLayout_InitialiseParticle;
   Particle_Index                                             totalInitialParticles = 0;
   double                                            averageInitialParticlesPerCell = 0.0;

   return (void*)_VoxelParticleLayout_New(  VOXELPARTICLELAYOUT_PASSARGS  );
}

void _VoxelParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data ) {
   VoxelParticleLayout*        self     = (VoxelParticleLayout*) particleLayout;
   Bool                        mustUseAllParticles;
   char*                       swarmVariableName;
   VoxelDataHandler_Abstract*  voxelDataHandler;
   Mesh*                       mesh;

   _GlobalParticleLayout_AssignFromXML( self, cf, data );
   /** Determines whether all particles must be placed within the domain.   If true, if any particles are not loaded, simulation will bail. */
   mustUseAllParticles = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"mustUseAllParticles", False );
   /** name of swarmvariable which data will be written to */
   swarmVariableName   = StG_Strdup( Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"swarmVariableName", "" ));
   /** voxel data handler to use */
   voxelDataHandler    = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"VoxelDataHandler", VoxelDataHandler_Abstract, True, data  );
   /** the mesh is an optional parameter so that we can make a good a priori estimate of required particles. */
   /** this is useful to avoid potentially prohibitive memory allocation attempts  */
   mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"GeometryMesh", Mesh, False, data  );

   Journal_Firewall( strcmp(swarmVariableName, ""),
        Journal_MyStream( Error_Type, self ),
        "\n\nError in %s for %s '%s' - No swarm variable name provided.\nIf you are using an rbf swarm, you will probably want to use something like 'mySwarm_rbfParticleSwarm-ParticleData' where \n'mySwarm' is what you have named your swarm in xml.\n\n\n",
        __func__,
        self->type,
        self->name );

   _VoxelParticleLayout_Init( self, mustUseAllParticles, swarmVariableName, voxelDataHandler, mesh );

}

void _VoxelParticleLayout_Init( void* particleLayout, Bool mustUseAllParticles, char* swarmVariableName, VoxelDataHandler_Abstract* voxelDataHandler, Mesh* geometryMesh )
{
   VoxelParticleLayout* self = (VoxelParticleLayout*) particleLayout;

   self->mustUseAllParticles = mustUseAllParticles;
   self->errorStream         = Journal_MyStream( Error_Type, self );
   self->swarmVariableName   = swarmVariableName;
   self->voxelDataHandler    = voxelDataHandler;
   self->geometryMesh        = geometryMesh;
}

void _VoxelParticleLayout_Build( void* particleLayout, void* data ) {
   VoxelParticleLayout* self = (VoxelParticleLayout*) particleLayout;
   Stg_Component_Build( self->voxelDataHandler, data, False );
   if( self->geometryMesh ) Stg_Component_Build( self->geometryMesh    , data, False );
}

void _VoxelParticleLayout_Initialise( void* particleLayout, void* data ) {
   VoxelParticleLayout* self = (VoxelParticleLayout*) particleLayout;
   Stg_Component_Initialise( self->voxelDataHandler, data, False );
   if( self->geometryMesh ) Stg_Component_Initialise( self->geometryMesh    , data, False );
}

void _VoxelParticleLayout_Execute( void* particleLayout, void* data ) {
   VoxelParticleLayout* self = (VoxelParticleLayout*) particleLayout;
   Stg_Component_Execute( self->voxelDataHandler, data, False );
   if( self->geometryMesh ) Stg_Component_Execute( self->geometryMesh    , data, False );
}

void _VoxelParticleLayout_Destroy( void* particleLayout, void* data ) {
   VoxelParticleLayout*  self = (VoxelParticleLayout*)particleLayout;

   Stg_Component_Destroy( self->voxelDataHandler, data, False );
   if( self->geometryMesh ) Stg_Component_Destroy( self->geometryMesh    , data, False );
   _GlobalParticleLayout_Destroy( self, data );
}

void _VoxelParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm ) {
   VoxelParticleLayout* self      = (VoxelParticleLayout*)particleLayout;
   Swarm*               swarm     = (Swarm*)_swarm;

   Journal_DPrintf( self->debug, "In %s(): for ParticleLayout \"%s\", of type %s\n",
      __func__, self->name, self->type );
   Stream_IndentBranch( Swarm_Debug );

   VoxelDataHandler_GetTotalVoxelCount( self->voxelDataHandler, &(self->totalInitialParticles) );

   _GlobalParticleLayout_SetInitialCounts( self, swarm );

   Stream_UnIndentBranch( Swarm_Debug );
   Journal_DPrintf( self->debug, "...finished %s() for ParticleLayout \"%s\".\n",
      __func__, self->name );
}

void _VoxelParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm ) {
   VoxelParticleLayout*                       self = (VoxelParticleLayout*)particleLayout;
   Swarm*                                    swarm = (Swarm*)_swarm;

   /** first get swarm variable which we will be initialising data on */
   self->swarmVar = SwarmVariable_Register_GetByName( swarm->swarmVariable_Register, self->swarmVariableName );
   Journal_Firewall( self->swarmVar,
        self->errorStream,
        "\n\nError in %s for %s '%s' - No swarm variable of the name '%s' found.  You may use the SwarmVariableList plugin to find which variables are available. \n\n\n",
        __func__,
        self->type,
        self->name,
        self->swarmVariableName );

   /** initialise voxel data handler */
   VoxelDataHandler_GotoVoxelDataTop( self->voxelDataHandler );
   
   _GlobalParticleLayout_InitialiseParticles( self, _swarm );


}

void _VoxelParticleLayout_InitialiseParticle( void* particleLayout, void* _swarm, Particle_Index newParticle_I, void* _particle ){
   VoxelParticleLayout*               self = (VoxelParticleLayout*)particleLayout;
   Swarm                            *swarm = (Swarm*)_swarm;
   GlobalParticle* particle = (GlobalParticle*)_particle;
   
   /** lets set the position */
   VoxelDataHandler_GetCurrentVoxelCoord( self->voxelDataHandler, &(particle->coord) );

   /** lets set the data */
   VoxelDataHandler_GetCurrentVoxelData(  self->voxelDataHandler, self->voxelData, self->swarmVar->variable->dataTypes[0] );
   Variable_SetValue( self->swarmVar->variable, swarm->particleLocalCount, self->voxelData );

   /** increment voxel data cursor */
   VoxelDataHandler_IncrementVoxelIndex( self->voxelDataHandler );
}

