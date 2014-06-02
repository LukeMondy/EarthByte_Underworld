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

#ifndef __ImportersToolbox_Base_VoxelFieldVariable_h__
#define __ImportersToolbox_Base_VoxelFieldVariable_h__

   /** Textual name of this class */
   extern const Type VoxelFieldVariable_Type;

   /* typedefs for function pointers */
   typedef void (VoxelFieldVariable_SetupDataArray)   ( void* voxelFieldVariable );
   typedef void (VoxelFieldVariable_TestDataArray)    ( void* voxelFieldVariable );
   typedef void (VoxelFieldVariable_SetArrayValue)    ( void* voxelFieldVariable, double* coord);
   typedef Bool (VoxelFieldVariable_GetArrayIndexIJK) ( void* voxelFieldVariable, Coord coord, int* indexIJK );

   /** VoxelFieldVariable contents */
   #define __VoxelFieldVariable                                              \
      __FieldVariable                                                        \
      VoxelFieldVariable_SetupDataArray*          _setupDataArray;          /** voxelDataArray configuration (determine size/type and build) routine */           \
      VoxelFieldVariable_TestDataArray*           _testDataArray;           /** auxiliary function to test data has been initialised correctly */                 \
      VoxelFieldVariable_SetArrayValue*           _setArrayValue;           /** retrieves values stored at provided voxelDataArray index */                       \
      VoxelFieldVariable_GetArrayIndexIJK*        _getArrayIndexIJK;        /** gets index into the voxelDataArray at a particular coordinate */                  \
      VoxelDataHandler_Abstract*  voxelDataHandler;                         /** data handler uses to handle voxel data */                                         \
      Mesh*                       decompositionMesh;                        /** mesh used to determine parallel decomposition (and so define only local data) */  \
      double                      localCrdMin[3];                           /** local coordinate minimums */                                                      \
      double                      localCrdMax[3];                           /** local coordinate maximums */                                                      \
      double                      localStartCrd[3];                         /** local coordinate starting position */                                             \
      unsigned                    arraySize[3];                             /** number of voxels in each direction */                                             \
      void*                       voxelDataArray;                           /** 3d array used to store data */                                                    \
      Variable_DataType           dataType;                                 /** data array dataType */                                                            \
      double                      localMin;                                 /** minimum value stored in the array locally */                                      \
      double                      localMax;                                 /** maximum value stored in the array locally */                                      \
      double                      zSliceCoord;                              /** z location of 3D voxel data to take slice from for 2d simulations */              \
      Bool                        useNearestCell;                           /** if querying outside the domain, should we used nearest cell? */                   \
      Stream*                     errorStream; \
      Bool                        minMaxCached;


   struct VoxelFieldVariable { __VoxelFieldVariable };

   /** Creation implementation */
   void* _VoxelFieldVariable_DefaultNew( Name name );


   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELFIELDVARIABLE_DEFARGS                                                        \
               FIELDVARIABLE_DEFARGS,                                                        \
                  VoxelFieldVariable_SetupDataArray*           _setupDataArray,              \
                  VoxelFieldVariable_TestDataArray*             _testDataArray,              \
                  VoxelFieldVariable_SetArrayValue*            _setArrayValue,               \
                  VoxelFieldVariable_GetArrayIndexIJK*         _getArrayIndexIJK


   #define VOXELFIELDVARIABLE_PASSARGS      \
               FIELDVARIABLE_PASSARGS,      \
                  _setupDataArray,          \
                  _testDataArray,           \
                  _setArrayValue,           \
                  _getArrayIndexIJK

   VoxelFieldVariable* _VoxelFieldVariable_New(  VOXELFIELDVARIABLE_DEFARGS  );

   /** Member initialisation implementation */

   void _VoxelFieldVariable_AssignFromXML( void* VoxelFieldVariable, Stg_ComponentFactory* cf, void* data ) ;
   void _VoxelFieldVariable_Build( void* VoxelFieldVariable, void* data ) ;
   #define _VoxelFieldVariable_Execute _FieldVariable_Execute
   void _VoxelFieldVariable_Destroy( void* VoxelFieldVariable, void* data ) ;
   void _VoxelFieldVariable_Initialise( void* VoxelFieldVariable, void* data ) ;

   InterpolationResult _VoxelFieldVariable_InterpolateValueAt( void* voxelFieldVariable, double* coord, double* value );
   double _VoxelFieldVariable_GetMinGlobalFieldMagnitude( void* voxelFieldVariable );
   double _VoxelFieldVariable_GetMaxGlobalFieldMagnitude( void* voxelFieldVariable );
   void _VoxelFieldVariable_CacheMinMaxGlobalFieldMagnitude( void* voxelFieldVariable );
   void _VoxelFieldVariable_GetMinAndMaxLocalCoords( void* voxelFieldVariable, double* min, double* max ) ;
   void _VoxelFieldVariable_GetMinAndMaxGlobalCoords( void* voxelFieldVariable, double* min, double* max ) ;


   void _VoxelFieldVariable_SetupDataArray( void* voxelFieldVariable );
   void _VoxelFieldVariable_TestDataArray( void* voxelFieldVariable );
   void _VoxelFieldVariable_SetArrayValue( void* self, double* coord );
   Bool _VoxelFieldVariable_GetArrayIndexIJK( void* self, Coord coord, int* indexIJK);
#endif /* __ImportersToolbox_Base_VoxelFieldVariable_h__ */

