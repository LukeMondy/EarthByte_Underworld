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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_GocadMaterials_h__
#define __ImportersToolbox_Base_VoxelDataHandler_GocadMaterials_h__

   /* Textual name of this class */
   extern const Type VoxelDataHandler_GocadMaterials_Type;

        /* VoxelDataHandler_GocadMaterials information */
   #define __VoxelDataHandler_GocadMaterials \
              __VoxelDataHandler_GocadAbstract \
                unsigned     numberRegions;         /** number of regions contained with gocad file */                       \
                unsigned     esize;                 /** actual stored size of datum */                                       \
                signed char* matIndexMapping;       /** maps from position in datum to a material Index */                   \
                HashTable*   hashTbl;               /** store hash table until materials have been instantiated */           \

   struct VoxelDataHandler_GocadMaterials { __VoxelDataHandler_GocadMaterials };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_GOCADMATERIALS_DEFARGS                                             \
      VOXELDATAHANDLER_GOCADABSTRACT_DEFARGS

   #define VOXELDATAHANDLER_GOCADMATERIALS_PASSARGS                                            \
      VOXELDATAHANDLER_GOCADABSTRACT_PASSARGS

   VoxelDataHandler_GocadMaterials* _VoxelDataHandler_GocadMaterials_New(  VOXELDATAHANDLER_GOCADMATERIALS_DEFARGS  );

   void _VoxelDataHandler_GocadMaterials_Init( VoxelDataHandler_GocadMaterials* voxelDataHandler );

   /* 'Stg_Class' Stuff */
   #define VoxelDataHandler_GocadMaterials_Copy( self ) \
      (VoxelDataHandler_GocadMaterials*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_GocadMaterials_DeepCopy( self ) \
      (VoxelDataHandler_GocadMaterials*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_GocadMaterials_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
   void _VoxelDataHandler_GocadMaterials_Delete( void* voxelDataHandler );

   /* 'Stg_Component' Stuff */
   void* _VoxelDataHandler_GocadMaterials_DefaultNew( Name name ) ;
   void _VoxelDataHandler_GocadMaterials_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_GocadMaterials_Build( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_GocadMaterials_Destroy( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_GocadMaterials_LoadMetaData( void* voxelDataHandler );
   int _VoxelDataHandler_GocadMaterials_GetCurrentVoxelData( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_GocadMaterials_h__ */

