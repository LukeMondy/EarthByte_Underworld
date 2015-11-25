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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_GocadProperties_h__
#define __ImportersToolbox_Base_VoxelDataHandler_GocadProperties_h__

   /* Textual name of this class */
   extern const Type VoxelDataHandler_GocadProperties_Type;

        /* VoxelDataHandler_GocadProperties information */
   #define __VoxelDataHandler_GocadProperties \
              __VoxelDataHandler_GocadAbstract \
                char      property_datatype[20]; /** datatype of property to be loaded */                                 \
                double    nodata_value;          /** value which signifies that no data exists */                         \
                char*     requiredPropertyName;  /** name of gocad property this handler will handle */                   \
                unsigned  requiredPropertyID;    /** alternatively the user can provide the ID of the required property */\

   struct VoxelDataHandler_GocadProperties { __VoxelDataHandler_GocadProperties };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_GOCADPROPERTIES_DEFARGS                                             \
      VOXELDATAHANDLER_GOCADABSTRACT_DEFARGS

   #define VOXELDATAHANDLER_GOCADPROPERTIES_PASSARGS                                            \
      VOXELDATAHANDLER_GOCADABSTRACT_PASSARGS

   VoxelDataHandler_GocadProperties* _VoxelDataHandler_GocadProperties_New(  VOXELDATAHANDLER_GOCADPROPERTIES_DEFARGS  );

   void _VoxelDataHandler_GocadProperties_Init( VoxelDataHandler_GocadProperties* voxelDataHandler, char* requiredPropertyName, unsigned requiredPropertyID );

   /* 'Stg_Class' Stuff */
   #define VoxelDataHandler_GocadProperties_Copy( self ) \
      (VoxelDataHandler_GocadProperties*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_GocadProperties_DeepCopy( self ) \
      (VoxelDataHandler_GocadProperties*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_GocadProperties_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void* _VoxelDataHandler_GocadProperties_DefaultNew( Name name ) ;
   void _VoxelDataHandler_GocadProperties_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_GocadProperties_Destroy( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_GocadProperties_LoadMetaData( void* voxelDataHandler );
   int _VoxelDataHandler_GocadProperties_GetCurrentVoxelData( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_GocadProperties_h__ */

