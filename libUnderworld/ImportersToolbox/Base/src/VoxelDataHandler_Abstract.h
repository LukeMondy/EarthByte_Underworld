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

#ifndef __ImportersToolbox_Base_VoxelDataHandler_Abstract_h__
#define __ImportersToolbox_Base_VoxelDataHandler_Abstract_h__

#define VOXELMAX_LINE_LENGTH_DEFINE 10000

   /* Textual name of this class */
   extern const Type VoxelDataHandler_Abstract_Type;

   /* typedefs for function pointers */
   typedef int (VoxelDataHandler_GetTotalVoxelCount)      ( void* voxelDataHandler, int* count );
   typedef int (VoxelDataHandler_GotoVoxelDataTop)        ( void* voxelDataHandler );
   typedef int (VoxelDataHandler_IncrementVoxelIndex)     ( void* voxelDataHandler );
   typedef int (VoxelDataHandler_GetCurrentVoxelIndex)    ( void* voxelDataHandler, int* index );
   typedef int (VoxelDataHandler_GetCurrentVoxelCoord)    ( void* voxelDataHandler, double* position );
   typedef int (VoxelDataHandler_GetCurrentVoxelDataFunc) ( void* voxelDataHandler );

   /* the NO_MATERIAL_TAG #def is used to indicate where no material should exist */
   #define NO_MATERIAL_TAG -1

   /* VoxelDataHandler_Abstract information */
   #define __VoxelDataHandler_Abstract \
      __Stg_Component                                                                                                     \
      VoxelDataHandler_GetTotalVoxelCount*        _getTotalVoxelCount;                                                    \
      VoxelDataHandler_GotoVoxelDataTop*            _gotoVoxelDataTop;                                                    \
      VoxelDataHandler_IncrementVoxelIndex*      _incrementVoxelIndex;                                                    \
      VoxelDataHandler_GetCurrentVoxelIndex*    _getCurrentVoxelIndex;                                                    \
      VoxelDataHandler_GetCurrentVoxelCoord*    _getCurrentVoxelCoord;                                                    \
      VoxelDataHandler_GetCurrentVoxelDataFunc*  _getCurrentVoxelData;                                                    \
      AbstractContext*        context;                                                                                    \
      Name                    filename;     /** voxel filename.. opening will be handled within implementation classes */ \
      File*                   fileMeta;     /** Manages the c file pointer */                                             \
      Stream*                 errorStream;  /** error stream for general error reporting */                               \
      Index                   numCells[3];  /** number of cells in each direction */                                      \
      double                  startCrd[3];  /** start coordinate of the first voxel centre */                             \
      double                  voxelSize[3]; /** voxel cell size in each direction */                                      \
      double                  scale[3];     /** factor to scale the voxet by (in each direction) */                       \
      int                     cursor;       /** 1d cursor to keep track of which voxel we are in. */                      \
      int                     cursorIJK[3]; /** cursor to keep track of which voxel we are in.. index 1 moves fastest */  \
      void*                   tempDataBuffer;      /** buffer to temporarily store voxel data. use for binary reads */    \
      size_t                  tempDataBufferSize;  /** size of temp buffer */                                             \
      void*                   dataBuffer;      /** buffer to store voxel data */                                          \
      size_t                  dataBufferSize;  /** size of buffer */                                                      \
      Variable_DataType       dataType;            /** datatype returned by GetCurrentVoxelData function */               \
      double                  minVoxelCoordsXYZ[3];  /** minimum coordinates for set of voxel centroids */                \
      double                  maxVoxelCoordsXYZ[3];  /** maximum coordinates for set of voxel centroids */                \
      double                  minDomainCoordsXYZ[3]; /** minimum coordinates for voxel domain */                          \
      double                  maxDomainCoordsXYZ[3]; /** maximum coordinates for voxel domain */                          \
      Index                   switchAxis[3];/** handles switching of axis */                                              \
      char                    currentString[VOXELMAX_LINE_LENGTH_DEFINE];  /** current string buffer */                   \
      char*                   currentToken;       /** current token */                                                    \
      char*                   currentStrtokKey;   /** key for strtok.  safer, as allows concurrent strtok usage */        \
      char                    delim[32];          /** string which specifies deliminators */                              \
      int                     dataStride;         /** stride size between datums */                                       \
      int                     dataPos;            /** position of required voxel data in data */                          \
      unsigned                bigEndian;          /** for binary data, signifies if data is bigendian */                  \
      unsigned                cursorData;         /** stores current 1D cursor position for data residing in the buffer */\
      void*                   currentData;        /** data from current cursor position for binary reading */             \
      File*                   fileData;           /** managed c file pointer to heavy data file */





   struct VoxelDataHandler_Abstract { __VoxelDataHandler_Abstract };

   /* Creation implementation / Virtual constructor */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define VOXELDATAHANDLER_ABSTRACT_DEFARGS                                             \
      STG_COMPONENT_DEFARGS,                                                             \
      VoxelDataHandler_GetTotalVoxelCount*                                   _getTotalVoxelCount,  \
      VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop,  \
      VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex,  \
      VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex,  \
      VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord,  \
      VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData

   #define VOXELDATAHANDLER_ABSTRACT_PASSARGS   \
      STG_COMPONENT_PASSARGS,                   \
         _getTotalVoxelCount,                        \
         _gotoVoxelDataTop,                     \
         _incrementVoxelIndex,                  \
         _getCurrentVoxelIndex,                 \
         _getCurrentVoxelCoord,                 \
         _getCurrentVoxelData

   VoxelDataHandler_Abstract* _VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_DEFARGS  );

   void _VoxelDataHandler_Abstract_Init( void* voxelDataHandler, AbstractContext* context, Name filename, char iAxisMap, char jAxisMap, char kAxisMap, double* scale );
   /* 'Stg_Class' Stuff */
   void _VoxelDataHandler_Abstract_Delete( void* voxelDataHandler );
   void _VoxelDataHandler_Abstract_Print( void* voxelDataHandler, Stream* stream );
   #define VoxelDataHandler_Abstract_Copy( self ) \
      (VoxelDataHandler_Abstract*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
   #define VoxelDataHandler_Abstract_DeepCopy( self ) \
      (VoxelDataHandler_Abstract*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
   void* _VoxelDataHandler_Abstract_Copy( void* voxelDataHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

   /* 'Stg_Component' Stuff */
   void _VoxelDataHandler_Abstract_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data );
   void _VoxelDataHandler_Abstract_Build( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_Abstract_Initialise( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_Abstract_Execute( void* voxelDataHandler, void* data );
   void _VoxelDataHandler_Abstract_Destroy( void* voxelDataHandler, void* data );

   /* --- Virtual Function Implementations --- */
   #define VoxelDataHandler_GetTotalVoxelCount( voxelDataHandler, count )      ( ((VoxelDataHandler_Abstract*) voxelDataHandler)->_getTotalVoxelCount( voxelDataHandler, count ) )
   #define VoxelDataHandler_GotoVoxelDataTop( voxelDataHandler )               ( ((VoxelDataHandler_Abstract*) voxelDataHandler)->_gotoVoxelDataTop( voxelDataHandler ) )
   #define VoxelDataHandler_IncrementVoxelIndex( voxelDataHandler )            ( ((VoxelDataHandler_Abstract*) voxelDataHandler)->_incrementVoxelIndex( voxelDataHandler ) )
   #define VoxelDataHandler_GetCurrentVoxelIndex( voxelDataHandler, index )    ( ((VoxelDataHandler_Abstract*) voxelDataHandler)->_getCurrentVoxelIndex( voxelDataHandler, index ) )
   #define VoxelDataHandler_GetCurrentVoxelCoord( voxelDataHandler, position ) ( ((VoxelDataHandler_Abstract*) voxelDataHandler)->_getCurrentVoxelCoord( voxelDataHandler, position ) )

   /* --- Wrapped Virtual Function Implementations --- */
   int VoxelDataHandler_GetCurrentVoxelData( void* voxelDataHandler, void* dataPointer, Variable_DataType dataType );
   /** minor routine to get prefix of string */
   char* _VoxelDataHandler_AbstractGetPrefix( char* filename );
   /** reset counters */
   int _VoxelDataHandler_Abstract_GotoVoxelDataTop( void* voxelDataHandler );

   /** returns the index of the current voxel the cursor is sitting on */
   int _VoxelDataHandler_Abstract_GetCurrentVoxelIndex( void* voxelDataHandler, int* index );

   /** returns the total count of voxels */
   int _VoxelDataHandler_Abstract_GetTotalVoxelCount( void* voxelDataHandler, int* count );

   /** returns the coordinate of the current voxel */
   int _VoxelDataHandler_Abstract_GetCurrentVoxelCoord( void* voxelDataHandler, double* coord );

   /** increments the current voxel UVW coord */
   int _VoxelDataHandler_Abstract_IncrementVoxelIndex( void* voxelDataHandler );

   /** compelete setup */
   void _VoxelDataHandler_Abstract_CompleteSetup( void* voxelDataHandler );

   /* routine to calculate domain min and max */
   void _VoxelDataHandler_Abstract_CalcMinMaxDomainCoords( void* voxelDataHandler );
   void _VoxelDataHandler_Abstract_CalcMinMaxVoxelCentroidCoords( void* voxelDataHandler );

   void _VoxelDataHandler_Abstract_GetASCIIVoxelData(  void* voxelDataHandler );
   void _VoxelDataHandler_Abstract_GetBinaryVoxelData( void* voxelDataHandler );

   void _VoxelDataHandler_Abstract_GetNextValidASCIIToken( void* voxelDataHandler );

#endif /* __ImportersToolbox_Base_VoxelDataHandler_Abstract_h__ */

