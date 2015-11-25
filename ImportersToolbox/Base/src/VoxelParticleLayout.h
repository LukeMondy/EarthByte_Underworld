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
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** TODO

* generalise ascii reading to allow different deliminators
* more robust handling of black spaces / lines
* read in header from file
* swarmVariables have have a multilevel tree of parents.  currently only check one above for checkpointing.
*
*/

#ifndef __ImportersToolbox_Base_VoxelParticleLayout_h__
#define __ImportersToolbox_Base_VoxelParticleLayout_h__


   /* Textual name of this class */
   extern const Type VoxelParticleLayout_Type;

   /* VoxelParticleLayout information */
#define __VoxelParticleLayout                          \
      __GlobalParticleLayout                           \
      Stream*                     errorStream;         \
      char*                       swarmVariableName;   \
      VoxelDataHandler_Abstract*  voxelDataHandler;    \
      SwarmVariable*              swarmVar;            \
      Mesh*                       geometryMesh;        \
      void*                       voxelData;

   struct VoxelParticleLayout { __VoxelParticleLayout };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELPARTICLELAYOUT_DEFARGS        \
                GLOBALPARTICLELAYOUT_DEFARGS
   #define VOXELPARTICLELAYOUT_PASSARGS        \
                GLOBALPARTICLELAYOUT_PASSARGS

   VoxelParticleLayout* _VoxelParticleLayout_New(  VOXELPARTICLELAYOUT_DEFARGS  );

   void _VoxelParticleLayout_Init( void* particleLayout, Bool mustUseAllParticles, char* swarmVariableName, VoxelDataHandler_Abstract* voxelDataHandler, Mesh* geometryMesh );

   /* 'Stg_Class' Stuff */
   void _VoxelParticleLayout_Delete( void* particleLayout );
   void _VoxelParticleLayout_Print( void* particleLayout, Stream* stream );
   #define VoxelParticleLayout_Copy( self ) \
      (VoxelParticleLayout*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelParticleLayout_DeepCopy( self ) \
      (VoxelParticleLayout*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelParticleLayout_Copy( void* particleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void* _VoxelParticleLayout_DefaultNew( Name name ) ;
   void _VoxelParticleLayout_AssignFromXML( void* particleLayout, Stg_ComponentFactory *cf, void* data );
   void _VoxelParticleLayout_Build( void* particleLayout, void* data );
   void _VoxelParticleLayout_Initialise( void* particleLayout, void* data );
   void _VoxelParticleLayout_Execute( void* particleLayout, void* data );
   void _VoxelParticleLayout_Destroy( void* particleLayout, void* data );

   void _VoxelParticleLayout_SetInitialCounts( void* particleLayout, void* _swarm ) ;
   void _VoxelParticleLayout_InitialiseParticles( void* particleLayout, void* _swarm ) ;
   void _VoxelParticleLayout_InitialiseParticle( void* particleLayout, void* swarm, Particle_Index newParticle_I, void* particle);

#endif /* __ImportersToolbox_Base_VoxelParticleLayout_h__ */

