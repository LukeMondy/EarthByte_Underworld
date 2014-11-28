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

#include <libgen.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_Abstract_Type = "VoxelDataHandler_Abstract";


VoxelDataHandler_Abstract* _VoxelDataHandler_Abstract_New( VOXELDATAHANDLER_ABSTRACT_DEFARGS )
{
   VoxelDataHandler_Abstract* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_Abstract ) );
   self = (VoxelDataHandler_Abstract*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );

   /* set default attributes */
   self->_getTotalVoxelCount    = _getTotalVoxelCount;
   self->_gotoVoxelDataTop      = _gotoVoxelDataTop;
   self->_incrementVoxelIndex   = _incrementVoxelIndex;
   self->_getCurrentVoxelIndex  = _getCurrentVoxelIndex;
   self->_getCurrentVoxelCoord  = _getCurrentVoxelCoord;
   self->_getCurrentVoxelData   = _getCurrentVoxelData;
   self->dataBuffer             = NULL;
   self->dataBufferSize         = -1;
   self->fileMeta               = NULL;
   self->fileData               = NULL;

   /** set counters to zero */
   self->cursor       = 0;
   self->cursorIJK[0] = 0;
   self->cursorIJK[1] = 0;
   self->cursorIJK[2] = 0;

   /** we default to use axis switch such that:
       I, the fastest moving coord = x
       J, the second fastest coord = z
       K,        the slowest coord = y                  */
   self->switchAxis[0] = 0;
   self->switchAxis[1] = 2;
   self->switchAxis[2] = 1;

   /** some things for ascii data.  these are not always required */
   self->currentToken     = NULL;
   self->currentStrtokKey = NULL;
   self->dataStride       = -1;
   self->dataPos          = 1;
   strcpy( self->delim, " ,\t\r\n");

   /** some things for binary data.  these are not always required */
   self->tempDataBuffer         = NULL;
   self->tempDataBufferSize     = -1;
   self->bigEndian = 0;   /** default to little endian */

   return self;
}

void _VoxelDataHandler_Abstract_CompleteSetup( void* voxelDataHandler ){
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*)voxelDataHandler;
   int ii;
   /** ensure voxelSize is always positive */
   for(ii=0; ii<3; ii++){
      if(self->voxelSize[ii]<0){
         self->scale[ii]     = -self->scale[ii];     /** reverse the scale so that we count in opposite direction */
         self->voxelSize[ii] = -self->voxelSize[ii]; /** set voxelSize to be positive */
      }
   }

   _VoxelDataHandler_Abstract_CalcMinMaxDomainCoords( self );
   _VoxelDataHandler_Abstract_CalcMinMaxVoxelCentroidCoords( self );
}

void _VoxelDataHandler_Abstract_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   if( self->dataBuffer ) {
      Memory_Free( self->dataBuffer );
      self->dataBuffer = NULL;
   }
   if( self->tempDataBuffer ) {
      Memory_Free( self->tempDataBuffer );
      self->tempDataBuffer = NULL;
   }
   if( self->filename ) {
      Memory_Free( self->filename );
      self->filename = NULL;
   }
   if( self->fileMeta ) {
      Stg_Class_Delete( self->fileMeta );
      self->fileMeta = NULL;
   }
   if( self->fileData ) {
      Stg_Class_Delete( self->fileData );
      self->fileData = NULL;
   }

   /* Stg_Class_Delete parent class */
   _Stg_Component_Delete( self );
}

void _VoxelDataHandler_Abstract_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_Abstract (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _Stg_Component_Print( self, stream );

   /* VoxelDataHandler_Abstract */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void _VoxelDataHandler_Abstract_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_Abstract*    self     = (VoxelDataHandler_Abstract*) voxelDataHandler;
   Name                    filename;
   AbstractContext*        context;
   double                  scale[3];
   char                    iAxisMap,jAxisMap,kAxisMap;

   context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", AbstractContext, False, data );
   if( !context  )
      context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, False, data  );

   /** The name of the file containing the voxel data */
   filename = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"filename", NULL ));
   /** add the file dirname to the search paths.. useful for when we open files listed in meta files */
   File_AddPath(dirname(StG_Strdup(filename)));
   /** organise axis mapping.. note that by default J-->Z  K-->Y */
   iAxisMap = *Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"mapIAxisToStgAxis", "X" );
   jAxisMap = *Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"mapJAxisToStgAxis", "Z" );
   kAxisMap = *Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"mapKAxisToStgAxis", "Y" );

   /** Scale IJK Axis  */
   scale[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"IAxisScale", 1 );
   scale[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"JAxisScale", 1 );
   scale[2] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"KAxisScale", 1 );

   _VoxelDataHandler_Abstract_Init( self, context, filename, iAxisMap, jAxisMap, kAxisMap, scale );

}

void _VoxelDataHandler_Abstract_Init( void* voxelDataHandler, AbstractContext* context, Name filename, char iAxisMap, char jAxisMap, char kAxisMap, double* scale )
{
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*) voxelDataHandler;

   self->context      = context;
   self->filename     = StG_Strdup( filename );
   self->errorStream  = Journal_MyStream( Error_Type, self );
   self->scale[0]     = scale[0];
   self->scale[1]     = scale[1];
   self->scale[2]     = scale[2];

   switch ( iAxisMap )
   {
      case 'X'  : self->switchAxis[0] = 0;  break;
      case 'Y'  : self->switchAxis[1] = 0;  break;
      case 'Z'  : self->switchAxis[2] = 0;  break;
      case 'x'  : self->switchAxis[0] = 0;  break;
      case 'y'  : self->switchAxis[1] = 0;  break;
      case 'z'  : self->switchAxis[2] = 0;  break;
      default   : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - 'iAxisMap' can only take values of 'X', 'Y' or 'Z'.",__func__,self->type,self->name ); break;
   }

   switch ( jAxisMap )
   {
      case 'X'  : self->switchAxis[0] = 1;  break;
      case 'Y'  : self->switchAxis[1] = 1;  break;
      case 'Z'  : self->switchAxis[2] = 1;  break;
      case 'x'  : self->switchAxis[0] = 1;  break;
      case 'y'  : self->switchAxis[1] = 1;  break;
      case 'z'  : self->switchAxis[2] = 1;  break;
      default   : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - 'jAxisMap' can only take values of 'X', 'Y' or 'Z'.",__func__,self->type,self->name ); break;
   }

   switch ( kAxisMap )
   {
      case 'X'  : self->switchAxis[0] = 2;  break;
      case 'Y'  : self->switchAxis[1] = 2;  break;
      case 'Z'  : self->switchAxis[2] = 2;  break;
      case 'x'  : self->switchAxis[0] = 2;  break;
      case 'y'  : self->switchAxis[1] = 2;  break;
      case 'z'  : self->switchAxis[2] = 2;  break;
      default   : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - 'kAxisMap' can only take values of 'X', 'Y' or 'Z'.",__func__,self->type,self->name ); break;
   }

   if( (self->switchAxis[0] == self->switchAxis[1]) || (self->switchAxis[0] == self->switchAxis[2]) || (self->switchAxis[1] == self->switchAxis[2]) )
      Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Axis mapping is not mutually exclusive! Values are %i, %i, %i.",__func__,self->type,self->name, self->switchAxis[0], self->switchAxis[1], self->switchAxis[2] );
}

void _VoxelDataHandler_Abstract_Build( void* voxelDataHandler, void* data ) {}

void _VoxelDataHandler_Abstract_Initialise( void* voxelDataHandler, void* data ) {}

void _VoxelDataHandler_Abstract_Execute( void* voxelDataHandler, void* data ) {}

void _VoxelDataHandler_Abstract_Destroy( void* voxelDataHandler, void* data ) {}

int _VoxelDataHandler_Abstract_GetCurrentVoxelIndex( void* voxelDataHandler, int* index ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   *index = self->cursor;
   return 1;
}

int _VoxelDataHandler_Abstract_GetTotalVoxelCount( void* voxelDataHandler, int* count ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   *count = self->numCells[0] * self->numCells[1] * self->numCells[2];
   return 1;
}

int _VoxelDataHandler_Abstract_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   /* init cursor */
   self->cursor       = 0;
   self->cursorData   = 0;
   self->cursorIJK[0] = 0;
   self->cursorIJK[1] = 0;
   self->cursorIJK[2] = 0;

   return 1;
}

int _VoxelDataHandler_Abstract_GetCurrentVoxelCoord( void* voxelDataHandler, double* coord ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;
   Coord voxelSystemCoord;

   /** get coordinates in voxel coord system */
   voxelSystemCoord[0] = self->startCrd[0] + self->cursorIJK[0]*self->voxelSize[0]*self->scale[0];
   voxelSystemCoord[1] = self->startCrd[1] + self->cursorIJK[1]*self->voxelSize[1]*self->scale[1];
   voxelSystemCoord[2] = self->startCrd[2] + self->cursorIJK[2]*self->voxelSize[2]*self->scale[2];

   /** now translate */
   coord[0] = voxelSystemCoord[self->switchAxis[0]];
   coord[1] = voxelSystemCoord[self->switchAxis[1]];
   coord[2] = voxelSystemCoord[self->switchAxis[2]];
   return 1;
}

int _VoxelDataHandler_Abstract_IncrementVoxelIndex( void* voxelDataHandler ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   self->cursor++;
   self->cursorIJK[0]++;
   /** if I cursor is past end, set it to zero, increment J, and read next line, set token cursor to zero */
   if(self->cursorIJK[0] == self->numCells[0]){
      self->cursorIJK[0] = 0;
      self->cursorIJK[1]++;
   }
   /** if J cursor is past end, set it to zero, increment K */
   if(self->cursorIJK[1] == self->numCells[1]){
      self->cursorIJK[1] = 0;
      self->cursorIJK[2]++;
   }
   return 1;
}

void _VoxelDataHandler_Abstract_CalcMinMaxVoxelCentroidCoords( void* voxelDataHandler ){
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*) voxelDataHandler;
   double displacedCoord;
   double minIJK[3], maxIJK[3];
   unsigned ii;

   /** calc min max IJK coords */
   for(ii=0; ii<3; ii++){
      /** calculate start coord + voxel displacement */
      displacedCoord = self->startCrd[ii] + self->scale[ii]*self->numCells[ii]*self->voxelSize[ii];
      /** mins */
      minIJK[ii] = (self->startCrd[ii] < displacedCoord) ? self->startCrd[ii]: displacedCoord;
      /** max */
      maxIJK[ii] = (self->startCrd[ii] > displacedCoord) ? self->startCrd[ii]: displacedCoord;
   }
   /** now translate to XYZ */
   self->minVoxelCoordsXYZ[0] = minIJK[self->switchAxis[0]];
   self->minVoxelCoordsXYZ[1] = minIJK[self->switchAxis[1]];
   self->minVoxelCoordsXYZ[2] = minIJK[self->switchAxis[2]];
   self->maxVoxelCoordsXYZ[0] = maxIJK[self->switchAxis[0]];
   self->maxVoxelCoordsXYZ[1] = maxIJK[self->switchAxis[1]];
   self->maxVoxelCoordsXYZ[2] = maxIJK[self->switchAxis[2]];
}

void _VoxelDataHandler_Abstract_CalcMinMaxDomainCoords( void* voxelDataHandler ){
   VoxelDataHandler_Abstract* self = (VoxelDataHandler_Abstract*) voxelDataHandler;
   double displacedCoord, startDomainCoord;
   double minIJK[3], maxIJK[3];
   unsigned ii;

   /** calc min max IJK coords */
   for(ii=0; ii<3; ii++){
      /** calculate start coord + voxel displacement */
      startDomainCoord = self->startCrd[ii] + self->scale[ii]*(                  -0.5)*self->voxelSize[ii];
      displacedCoord   = self->startCrd[ii] + self->scale[ii]*(self->numCells[ii]-0.5)*self->voxelSize[ii];
      /** mins */
      minIJK[ii] = (startDomainCoord < displacedCoord) ? startDomainCoord: displacedCoord;
      /** max */
      maxIJK[ii] = (startDomainCoord > displacedCoord) ? startDomainCoord: displacedCoord;
   }
   /** now translate to XYZ */
   self->minDomainCoordsXYZ[0] = minIJK[self->switchAxis[0]];
   self->minDomainCoordsXYZ[1] = minIJK[self->switchAxis[1]];
   self->minDomainCoordsXYZ[2] = minIJK[self->switchAxis[2]];
   self->maxDomainCoordsXYZ[0] = maxIJK[self->switchAxis[0]];
   self->maxDomainCoordsXYZ[1] = maxIJK[self->switchAxis[1]];
   self->maxDomainCoordsXYZ[2] = maxIJK[self->switchAxis[2]];
}

int VoxelDataHandler_GetCurrentVoxelData( void* voxelDataHandler, void* dataPointer, Variable_DataType dataType ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;
   unsigned returnVal;
   /** call getvoxeldata function, which will set the value stored at self->dataBuffer */
   returnVal = self->_getCurrentVoxelData( self );

   /** copy to provided location, converting type if required */
   if( (self->dataType==dataType) || (dataType==(int)NULL) )
      memcpy(dataPointer, self->dataBuffer, self->dataBufferSize);

   else
      switch( dataType )
      {

      case Variable_DataType_Char :
         {
         signed char tempConvertedDataChar;
         switch ( self->dataType )
         {
            case Variable_DataType_Int    : tempConvertedDataChar = (signed char)(    *(int*)self->dataBuffer       ); break;
            case Variable_DataType_Float  : tempConvertedDataChar = (signed char)(  *(float*)self->dataBuffer + 0.5 ); break;
            case Variable_DataType_Double : tempConvertedDataChar = (signed char)( *(double*)self->dataBuffer + 0.5 ); break;

         }
         memcpy(dataPointer, (void*)&tempConvertedDataChar, sizeof(signed char));
         }
         break;

      case Variable_DataType_Int :
         {
         int tempConvertedDataInt;
         switch ( self->dataType )
         {
            case Variable_DataType_Char   : tempConvertedDataInt = (int)( *(signed char*)self->dataBuffer       ); break;
            case Variable_DataType_Float  : tempConvertedDataInt = (int)(       *(float*)self->dataBuffer + 0.5 ); break;
            case Variable_DataType_Double : tempConvertedDataInt = (int)(      *(double*)self->dataBuffer + 0.5 ); break;

         }
         memcpy(dataPointer, (void*)&tempConvertedDataInt, sizeof(int));
         }
         break;

      case Variable_DataType_Float :
         {
         float tempConvertedDataFloat;
         switch ( self->dataType )
         {
            case Variable_DataType_Char   : tempConvertedDataFloat = (float)( *(signed char*)self->dataBuffer ); break;
            case Variable_DataType_Int    : tempConvertedDataFloat = (float)(         *(int*)self->dataBuffer ); break;
            case Variable_DataType_Double : tempConvertedDataFloat = (float)(      *(double*)self->dataBuffer ); break;

         }
         memcpy(dataPointer, (void*)&tempConvertedDataFloat, sizeof(float));
         }
         break;

      case Variable_DataType_Double :
         {
         double tempConvertedDataDouble;
         switch ( self->dataType )
         {
            case Variable_DataType_Char  : tempConvertedDataDouble = (double)( *(signed char*)self->dataBuffer ); break;
            case Variable_DataType_Int   : tempConvertedDataDouble = (double)(         *(int*)self->dataBuffer ); break;
            case Variable_DataType_Float : tempConvertedDataDouble = (double)(       *(float*)self->dataBuffer ); break;

         }
         memcpy(dataPointer, (void*)&tempConvertedDataDouble, sizeof(double));
         }
         break;

      default :
         Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Datatype not know or not supported.",__func__,self->type,self->name ); break;

      }

      return returnVal;

}

void _VoxelDataHandler_Abstract_GetASCIIVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   while( self->cursorData <= self->cursor ){
      unsigned jj;  /* data position indicator */

      if ( self->dataStride > 0 ){                            /* for this mode (dataStride>0), we grab a new line per datum */
      /** in this case, each datum is on a new line, so keep grabbing new lines until we are up to current cursor position */
         char* result = fgets( &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
         if( self->currentString != result )
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s'.\n"
                                                    "Please ensure voxel dimensions are correct for this file and that the file is not corrupted.",__func__,self->type,self->name, self->fileData->name );
         /** if we are at the required position, tokenise */
         if( self->cursorData == self->cursor ){
            self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
            /** read to required position */
            for( jj=1; jj< self->dataPos; jj++ )
               self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
         }
      } else if ( self->dataStride < 0 ) {                     /* for this mode (dataStride<0), we simply keep reading, taking a new line as required */
         /** if this is our first entry since rewind, set token to null, skip forward to dataPos, else step forward by the stride length */
         if( self->cursor == 0 ) self->currentToken = NULL;
         int jjmax = self->cursor == 0 ? self->dataPos : abs(self->dataStride);
         for( jj = 0 ; jj<jjmax ; jj++)
            _VoxelDataHandler_Abstract_GetNextValidASCIIToken( voxelDataHandler );
      } else {
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - A datastride of zero is not valid. Please check your input.",__func__,self->type,self->name, self->fileData->name );
      }

      self->cursorData++;
   }

}

void _VoxelDataHandler_Abstract_GetNextValidASCIIToken( void* voxelDataHandler ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;

   if( self->currentToken == NULL ){
      while( self->currentToken == NULL ) {
         char* result = fgets( &self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
         if( self->currentString != result )
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s'.\n"
                                                    "Please ensure voxel dimensions are correct for this file and that the file is not corrupted.",__func__,self->type,self->name, self->fileData->name );
         self->currentToken = strtok_r( self->currentString, self->delim, &self->currentStrtokKey );
      }
   } else {
      self->currentToken = strtok_r( NULL, self->delim, &self->currentStrtokKey );
      if(self->currentToken == NULL)
         _VoxelDataHandler_Abstract_GetNextValidASCIIToken( voxelDataHandler );
   }

}

void _VoxelDataHandler_Abstract_GetBinaryVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_Abstract*  self = (VoxelDataHandler_Abstract*)voxelDataHandler;
   size_t resultSeek;
   size_t resultRead;

   if( self->cursorData <= self->cursor ){  /**need to first ensure file cursor is at correct position */
      char  turnaround[self->tempDataBufferSize];
      unsigned ii;

      resultSeek = fseek( CFile_Ptr( self->fileData ), (self->cursor - self->cursorData)*self->tempDataBufferSize, SEEK_CUR );
      self->cursorData = self->cursor;

      resultRead = fread( (void*)&turnaround, self->tempDataBufferSize, 1, CFile_Ptr( self->fileData ) );
      if( (resultRead != 1) || (resultSeek != 0) )
         Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s'.\n"
                                                 "Please ensure voxel dimensions are correct for this file and that the file is not corrupted.",__func__,self->type,self->name, self->fileData->name );

      if(self->bigEndian){
         /** covert to little-endian */
         for(ii=0; ii<self->tempDataBufferSize; ii++)
            ((char*)self->tempDataBuffer)[ii] = ((char*)&turnaround)[self->tempDataBufferSize-1-ii];
      } else {
         memcpy( self->tempDataBuffer, &turnaround ,self->tempDataBufferSize );
      }
      self->cursorData++;
   }
}

