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

#ifndef __ImportersToolbox_Base_VoxelFieldVariable_GMT_h__
#define __ImportersToolbox_Base_VoxelFieldVariable_GMT_h__

   /** Textual name of this class */
   extern const Type VoxelFieldVariable_GMT_Type;

   /** VoxelFieldVariable_GMT contents */
   #define __VoxelFieldVariable_GMT    \
           __VoxelFieldVariable        \
           Bool         cartesianMode; \
           Bool         mapTo360;      \
           double       plateDepth;    \
           double       innerRadius;


   struct VoxelFieldVariable_GMT { __VoxelFieldVariable_GMT };

   /** Creation implementation */
   void* _VoxelFieldVariable_GMT_DefaultNew( Name name );

   #define VOXELFIELDVARIABLE_GMT_DEFARGS        \
               VOXELFIELDVARIABLE_DEFARGS

   #define VOXELFIELDVARIABLE_GMT_PASSARGS      \
               VOXELFIELDVARIABLE_PASSARGS

   VoxelFieldVariable_GMT* _VoxelFieldVariable_GMT_New(  VOXELFIELDVARIABLE_GMT_DEFARGS  );

   /** Member initialisation implementation */

   void _VoxelFieldVariable_GMT_AssignFromXML( void* VoxelFieldVariable_GMT, Stg_ComponentFactory* cf, void* data ) ;
//Bool _VoxelFieldVariable_GMT_GetArrayIndexIJK( void* voxelFieldVariable, Coord coord, int* indexIJK);
   InterpolationResult _VoxelFieldVariable_GMT_InterpolateValueAt( void* voxelFieldVariable, double* coord, double* value );

#endif /* __ImportersToolbox_Base_VoxelFieldVariable_GMT_h__ */

