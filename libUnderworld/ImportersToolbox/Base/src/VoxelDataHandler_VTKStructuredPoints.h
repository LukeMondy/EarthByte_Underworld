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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_VTKStructuredPoints_h__
#define __ImportersToolbox_Base_VoxelDataHandler_VTKStructuredPoints_h__

#define STANDARDVTKFILECOMMENT "\nThe expected file should look very much like this:\n\n \
# vtk DataFile Version 2.0 \n \
Really cool data \n \
ASCII \n \
DATASET STRUCTURED_POINTS \n \
DIMENSIONS 3 3 3 \n \
ORIGIN -4 -4 -4 \n \
SPACING 2.66666666667 2.66666666667 2.66666666667 \n \
POINT_DATA 27 \n \
SCALARS points double 1 \n \
LOOKUP_TABLE default\n \
(data here)\n"

   /* Textual name of this class */
   extern const Type VoxelDataHandler_VTKStructuredPoints_Type;

        /* VoxelDataHandler_VTKStructuredPoints information */
   #define __VoxelDataHandler_VTKStructuredPoints \
   __VoxelDataHandler_Abstract \
      char*  datasetName;     \
      int    unsignedType;

   struct VoxelDataHandler_VTKStructuredPoints { __VoxelDataHandler_VTKStructuredPoints };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_VTKSTRUCTUREDPOINTS_DEFARGS                                             \
      VOXELDATAHANDLER_ABSTRACT_DEFARGS

   #define VOXELDATAHANDLER_VTKSTRUCTUREDPOINTS_PASSARGS                                            \
      VOXELDATAHANDLER_ABSTRACT_PASSARGS

   VoxelDataHandler_VTKStructuredPoints* _VoxelDataHandler_VTKStructuredPoints_New(  VOXELDATAHANDLER_VTKSTRUCTUREDPOINTS_DEFARGS  );

   /* 'Stg_Class' Stuff */
   void _VoxelDataHandler_VTKStructuredPoints_Delete( void* voxelDataHandler );
   void _VoxelDataHandler_VTKStructuredPoints_Print( void* voxelDataHandler, Stream* stream );
   #define VoxelDataHandler_VTKStructuredPoints_Copy( self ) \
      (VoxelDataHandler_VTKStructuredPoints*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_VTKStructuredPoints_DeepCopy( self ) \
      (VoxelDataHandler_VTKStructuredPoints*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_VTKStructuredPoints_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void* _VoxelDataHandler_VTKStructuredPoints_DefaultNew( Name name ) ;
   void _VoxelDataHandler_VTKStructuredPoints_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_VTKStructuredPoints_Build( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_VTKStructuredPoints_Initialise( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_VTKStructuredPoints_Execute( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_VTKStructuredPoints_Destroy( void* voxelDataHandler, void* data );

   int _VoxelDataHandler_VTKStructuredPoints_LoadMetaData( void* voxelDataHandler );
   int _VoxelDataHandler_VTKStructuredPoints_GotoVoxelDataTop( void* voxelDataHandler );

   void _VoxelDataHandler_VTKStructuredPoints_SkipBlank( File* file );
   int _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataASCII( void* voxelDataHandler );
   int _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataBinary( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_VTKStructuredPoints_h__ */

