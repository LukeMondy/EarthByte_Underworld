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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_ndarray_h__
#define __ImportersToolbox_Base_VoxelDataHandler_ndarray_h__

#include "VoxelDataHandler_ASCII.h"
   /* Textual name of this class */
   extern const Type VoxelDataHandler_ndarray_Type;

        /* VoxelDataHandler_ndarray information */
   #define __VoxelDataHandler_ndarray \
       __VoxelDataHandler_ASCII \
      size_t*                  ndPointer;         /** pointer to underlying numpy array memory */                  

   struct VoxelDataHandler_ndarray { __VoxelDataHandler_ndarray };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_NDARRAY_DEFARGS                                             \
      VOXELDATAHANDLER_ASCII_DEFARGS

   #define VOXELDATAHANDLER_NDARRAY_PASSARGS                                            \
      VOXELDATAHANDLER_ASCII_PASSARGS

   VoxelDataHandler_ndarray* _VoxelDataHandler_ndarray_New(  VOXELDATAHANDLER_NDARRAY_DEFARGS  );
   void _VoxelDataHandler_ndarray_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );

   /* 'Stg_Class' Stuff */
   #define VoxelDataHandler_ndarray_Copy( self ) \
      (VoxelDataHandler_ndarray*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_ndarray_DeepCopy( self ) \
      (VoxelDataHandler_ndarray*)Stg_Class_Copy( self, NULL, True, NULL, NULL )

   /* 'Stg_Component' Stuff */
   void* _VoxelDataHandler_ndarray_DefaultNew( Name name ) ;
   void _VoxelDataHandler_ndarray_Build( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_ndarray_GetCurrentVoxelData( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_ndarray_h__ */

