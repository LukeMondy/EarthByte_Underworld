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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "VoxelDataHandler_Abstract.h"
#include "VoxelDataHandler_ASCII.h"
#include "VoxelDataHandler_ndarray.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_ndarray_Type = "VoxelDataHandler_ndarray";

VoxelDataHandler_ndarray* _VoxelDataHandler_ndarray_New( VOXELDATAHANDLER_NDARRAY_DEFARGS )
{
   VoxelDataHandler_ndarray* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_ndarray ) );
   self = (VoxelDataHandler_ndarray*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->dataType = Variable_DataType_Double;

   return self;
}

void* _VoxelDataHandler_ndarray_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_ndarray);
   Type                                                                        type = VoxelDataHandler_ndarray_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_Abstract_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_Abstract_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_ndarray_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_ndarray_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_ndarray_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_Abstract_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_Abstract_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_Abstract_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_Abstract_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_ndarray_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_ndarray_New( VOXELDATAHANDLER_NDARRAY_PASSARGS  );
}

void _VoxelDataHandler_ndarray_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_ndarray*    self     = (VoxelDataHandler_ndarray*) voxelDataHandler;
   char* dataType;

   /** config parent */
   _VoxelDataHandler_ASCII_AssignFromXML( voxelDataHandler, cf, data );

   /** see header file for info, if not self evident */
   sscanf(Stg_ComponentFactory_GetString(    cf, self->name, (Dictionary_Entry_Key)"ndPointer", "" ),"%p",(void **)&self->ndPointer) ;
}

void _VoxelDataHandler_ndarray_Build( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_ndarray*    self     = (VoxelDataHandler_ndarray*) voxelDataHandler;

   _VoxelDataHandler_Abstract_Build( voxelDataHandler, data);

   /** setup voxel data ready to read, incase gototop routine isn't called */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( self );

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : self->dataBufferSize = sizeof(signed char);  break;
      case Variable_DataType_Int    : self->dataBufferSize = sizeof(        int);  break;
      case Variable_DataType_Float  : self->dataBufferSize = sizeof(      float);  break;
      case Variable_DataType_Double : self->dataBufferSize = sizeof(     double);  break;
      default                       : Journal_Firewall( NULL, self->errorStream, "Error in %s for %s '%s' - Datatype not know or not supported.",__func__,self->type,self->name ); break;
   }
   self->dataBuffer = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "voxel data handler type" );
}

int _VoxelDataHandler_ndarray_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_ndarray*  self = (VoxelDataHandler_ndarray*)voxelDataHandler;
   double doubleData;

   /* copy from buffer into dataBuffer */
   memcpy(self->dataBuffer, self->ndPointer+self->dataBufferSize*self->cursorIJK[0]*self->cursorIJK[1]*self->cursorIJK[2], self->dataBufferSize);

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : doubleData = (double) *(signed char*)self->dataBuffer; break;
      case Variable_DataType_Int    : doubleData = (double) *(        int*)self->dataBuffer; break;
      case Variable_DataType_Float  : doubleData = (double) *(      float*)self->dataBuffer; break;
      case Variable_DataType_Double : doubleData = (double) *(     double*)self->dataBuffer; break;
   }
   if(doubleData > self->noDataValue)
      return 1;
   else
      return 0;

}

