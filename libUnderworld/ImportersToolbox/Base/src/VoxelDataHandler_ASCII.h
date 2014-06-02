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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_ASCII_h__
#define __ImportersToolbox_Base_VoxelDataHandler_ASCII_h__

   /* Textual name of this class */
   extern const Type VoxelDataHandler_ASCII_Type;

        /* VoxelDataHandler_ASCII information */
   #define __VoxelDataHandler_ASCII \
   __VoxelDataHandler_Abstract \
      Bool                    metadataRead;        /** signifies whether metadata has been read */                        \
      double                  noDataValue;         /** string which signifies voxels devoid of data */                    \

   struct VoxelDataHandler_ASCII { __VoxelDataHandler_ASCII };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_ASCII_DEFARGS                                             \
      VOXELDATAHANDLER_ABSTRACT_DEFARGS

   #define VOXELDATAHANDLER_ASCII_PASSARGS                                            \
      VOXELDATAHANDLER_ABSTRACT_PASSARGS

   VoxelDataHandler_ASCII* _VoxelDataHandler_ASCII_New(  VOXELDATAHANDLER_ASCII_DEFARGS  );

   /* 'Stg_Class' Stuff */
   void _VoxelDataHandler_ASCII_Delete( void* voxelDataHandler );
   void _VoxelDataHandler_ASCII_Print( void* voxelDataHandler, Stream* stream );
   #define VoxelDataHandler_ASCII_Copy( self ) \
      (VoxelDataHandler_ASCII*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_ASCII_DeepCopy( self ) \
      (VoxelDataHandler_ASCII*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_ASCII_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void* _VoxelDataHandler_ASCII_DefaultNew( Name name ) ;
   void _VoxelDataHandler_ASCII_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_ASCII_Build( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_ASCII_Initialise( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_ASCII_Execute( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_ASCII_Destroy( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_ASCII_LoadMetaData( void* voxelDataHandler );
   int _VoxelDataHandler_ASCII_GetTotalVoxelCount( void* voxelDataHandler, int* count );
   int _VoxelDataHandler_ASCII_GotoVoxelDataTop( void* voxelDataHandler );
   int _VoxelDataHandler_ASCII_GetCurrentVoxelIndex( void* voxelDataHandler, int* index );
   int _VoxelDataHandler_ASCII_GetCurrentVoxelCoord( void* voxelDataHandler, double* coord );
   int _VoxelDataHandler_ASCII_GetCurrentVoxelData( void* voxelDataHandler );

   void _VoxelDataHandler_ASCII_SkipBlank( File* file );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_ASCII_h__ */

