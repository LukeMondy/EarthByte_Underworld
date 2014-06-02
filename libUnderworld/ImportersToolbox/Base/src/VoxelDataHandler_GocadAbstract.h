/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**      Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**      Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**      Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**      Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**      Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**      Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_GocadAbstract_h__
#define __ImportersToolbox_Base_VoxelDataHandler_GocadAbstract_h__

   /* Textual name of this class */
   extern const Type VoxelDataHandler_GocadAbstract_Type;

        /* VoxelDataHandler_GocadAbstract information */
   #define __VoxelDataHandler_GocadAbstract \
              __VoxelDataHandler_Abstract \
                double    coord_axis[3][3];      /** vector determining direction of voxel dataset axis */                \
                double    axis_min[3];           /** position in vox coord system where data starts */                    \
                double    axis_max[3];           /** position in vox coord system where data ends */                      \
                Bool      zPositiveUp;           /** signifies if z axis is positive upwards */                           \
                char      filenameData[300];     /** filename for file storing heavy/binary voxel data  */    


   struct VoxelDataHandler_GocadAbstract { __VoxelDataHandler_GocadAbstract };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_GOCADABSTRACT_DEFARGS                                             \
      VOXELDATAHANDLER_ABSTRACT_DEFARGS

   #define VOXELDATAHANDLER_GOCADABSTRACT_PASSARGS                                            \
      VOXELDATAHANDLER_ABSTRACT_PASSARGS

   VoxelDataHandler_GocadAbstract* _VoxelDataHandler_GocadAbstract_New(  VOXELDATAHANDLER_GOCADABSTRACT_DEFARGS  );

   void _VoxelDataHandler_GocadAbstract_Init( VoxelDataHandler_GocadAbstract* voxelDataHandler );

   /* 'Stg_Class' Stuff */
   void _VoxelDataHandler_GocadAbstract_Delete( void* voxelDataHandler );
   void _VoxelDataHandler_GocadAbstract_Print( void* voxelDataHandler, Stream* stream );
   #define VoxelDataHandler_GocadAbstract_Copy( self ) \
      (VoxelDataHandler_GocadAbstract*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_GocadAbstract_DeepCopy( self ) \
      (VoxelDataHandler_GocadAbstract*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_GocadAbstract_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void _VoxelDataHandler_GocadAbstract_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_GocadAbstract_Build( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_GocadAbstract_Initialise( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_GocadAbstract_Execute( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_GocadAbstract_Destroy( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_GocadAbstract_LoadMetaData( void* voxelDataHandler );
   int _VoxelDataHandler_GocadAbstract_GotoVoxelDataTop( void* voxelDataHandler );
   int _VoxelDataHandler_GocadAbstract_IncrementVoxelIndex( void* voxelDataHandler );
   int _VoxelDataHandler_GocadAbstract_GetCurrentVoxelData( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_GocadAbstract_h__ */

