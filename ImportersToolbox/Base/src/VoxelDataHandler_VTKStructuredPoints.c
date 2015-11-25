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
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "VoxelDataHandler_Abstract.h"
#include "VoxelDataHandler_VTKStructuredPoints.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_VTKStructuredPoints_Type = "VoxelDataHandler_VTKStructuredPoints";

VoxelDataHandler_VTKStructuredPoints* _VoxelDataHandler_VTKStructuredPoints_New( VOXELDATAHANDLER_VTKSTRUCTUREDPOINTS_DEFARGS )
{
   VoxelDataHandler_VTKStructuredPoints* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_VTKStructuredPoints ) );
   self = (VoxelDataHandler_VTKStructuredPoints*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->dataStride   = -1;
   self->dataPos      = 1;
   self->bigEndian    = 0;
   self->datasetName  = NULL;
   self->unsignedType = 0;

   return self;
}

void _VoxelDataHandler_VTKStructuredPoints_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_VTKStructuredPoints* self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;

   /* delete parent class */
   _VoxelDataHandler_Abstract_Delete( self );
}

void _VoxelDataHandler_VTKStructuredPoints_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_VTKStructuredPoints* self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_VTKStructuredPoints (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _VoxelDataHandler_Abstract_Print( self, stream );

   /* VoxelDataHandler_VTKStructuredPoints */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void* _VoxelDataHandler_VTKStructuredPoints_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_VTKStructuredPoints);
   Type                                                                        type = VoxelDataHandler_VTKStructuredPoints_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_VTKStructuredPoints_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_VTKStructuredPoints_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_VTKStructuredPoints_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_VTKStructuredPoints_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_VTKStructuredPoints_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_VTKStructuredPoints_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_VTKStructuredPoints_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_VTKStructuredPoints_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_VTKStructuredPoints_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = NULL;  /** set below */

   return (void*)_VoxelDataHandler_VTKStructuredPoints_New( VOXELDATAHANDLER_VTKSTRUCTUREDPOINTS_PASSARGS  );
}

void _VoxelDataHandler_VTKStructuredPoints_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_VTKStructuredPoints*    self     = (VoxelDataHandler_VTKStructuredPoints*) voxelDataHandler;

   /** config parent */
   _VoxelDataHandler_Abstract_AssignFromXML( voxelDataHandler, cf, data );

   self->datasetName = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"DatasetName", NULL ));

   /** open file */
   self->fileMeta  = CFile_NewRead( self->filename );
   Journal_Firewall( self->fileMeta != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filename );
   /** open a second copy */
   self->fileData  = CFile_NewRead( self->filename );
   /** go ahead and load meta data */
   _VoxelDataHandler_VTKStructuredPoints_LoadMetaData( voxelDataHandler );

   File_Close( self->fileMeta );

   /** now that we've loaded required information, setup auxilliary quantities */
   _VoxelDataHandler_Abstract_CompleteSetup( voxelDataHandler );

   /** setup voxel data ready to read, incase gototop routine isn't called */
   _VoxelDataHandler_VTKStructuredPoints_GotoVoxelDataTop( self );

}

void _VoxelDataHandler_VTKStructuredPoints_Build( void* voxelDataHandler, void* data ) { 
   _VoxelDataHandler_Abstract_Build( voxelDataHandler, data); 
   VoxelDataHandler_VTKStructuredPoints*    self     = (VoxelDataHandler_VTKStructuredPoints*) voxelDataHandler;

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : self->tempDataBufferSize = sizeof(signed char);  break;
      case Variable_DataType_Int    : self->tempDataBufferSize = sizeof(        int);  break;
      case Variable_DataType_Float  : self->tempDataBufferSize = sizeof(      float);  break;
      case Variable_DataType_Double : self->tempDataBufferSize = sizeof(     double);  break;
      default                       : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Datatype not know or not supported.",__func__,self->type,self->name ); break;
   }
   self->tempDataBuffer = Memory_Alloc_Bytes_Unnamed( self->tempDataBufferSize , "voxel data handler type" );
   /** final buffer and temp buffer are same size */
   self->dataBufferSize = self->tempDataBufferSize;
   self->dataBuffer = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "voxel data handler type" );

}

void _VoxelDataHandler_VTKStructuredPoints_Initialise( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Initialise( voxelDataHandler, data); }

void _VoxelDataHandler_VTKStructuredPoints_Execute( void* voxelDataHandler, void* data )  { _VoxelDataHandler_Abstract_Execute( voxelDataHandler, data); }

void _VoxelDataHandler_VTKStructuredPoints_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_VTKStructuredPoints*  self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_Abstract_Destroy( self, data );
}

int _VoxelDataHandler_VTKStructuredPoints_LoadMetaData( void* voxelDataHandler ){
   VoxelDataHandler_VTKStructuredPoints*  self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;
   char* identifier;
   char* result;
   unsigned jj;
   int pointdata;
   /** rewind file to top */
   rewind( CFile_Ptr( self->fileMeta ) );
   /** skip blank lines (if any) */
   _VoxelDataHandler_VTKStructuredPoints_SkipBlank( self->fileMeta );

   for(jj=1; jj<=10 ;jj++){

      /** now read in data from file. go straight to third line */
      result = fgets( self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileMeta ) );
      if( self->currentString != result )
         Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s' before voxel meta information could be read.\n"
                                                 "Please ensure that the file is not corrupted.",__func__,self->type,self->name, self->fileMeta->name );
      identifier = strtok_r(self->currentString, self->delim, &self->currentStrtokKey );

      if(jj==3){
         if(      !strcmp(identifier, "ASCII"  )){
            self->_getCurrentVoxelData = _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataASCII;
         }
         else if (!strcmp(identifier, "BINARY" )){
            self->_getCurrentVoxelData = _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataBinary;
         }
         else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 3, expected key 'BINARY' or 'ASCII', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==4){
         if (!strcmp(identifier, "DATASET")){
            char* token = strtok_r(NULL, self->delim, &self->currentStrtokKey );
            if ( strcmp( token, "STRUCTURED_POINTS" ) )
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nVTK format '%s' is not supported.\n"
                                                       "Please ensure that your file is of VTK STRUCTURED_POINTS format.",__func__,self->type,self->name, token );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 4, expected key 'DATASET', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==5){
         if (!strcmp(identifier, "DIMENSIONS")){
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u", &self->numCells[0] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u", &self->numCells[1] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u", &self->numCells[2] );
            if( (self->numCells[0]<=0) || (self->numCells[1]<=0) || (self->numCells[2]<=0) )
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nVTK file dimensions of (%u,%u,%u) no valid.\n"
                                                       "Voxel dimensions should be >0",__func__,self->type,self->name, self->numCells[0],self->numCells[1],self->numCells[2] );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 5, expected key 'DIMENSIONS', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==6){
         if (!strcmp(identifier, "ORIGIN")){
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[0] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[1] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[2] );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 6, expected key 'ORIGIN', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==7){
         if (!strcmp(identifier, "SPACING") || !strcmp(identifier, "ASPECT_RATIO")  ){
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[0] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[1] );
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[2] );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 7, expected key 'ORIGIN' or 'ASPECT_RATIO', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==8){
         if (!strcmp(identifier, "POINT_DATA")){
            sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u", &pointdata );
            if( pointdata != self->numCells[0]*self->numCells[1]*self->numCells[2])
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nSeemingly inconsistent meta data encountered in file '%s'\n"
                                                       "POINT_DATA attribute suggests %u points, but file DIMENSIONS (%u,%u,%u) suggests %u points.\n ", \
                                                       __func__,self->type,self->name, self->fileMeta->name, pointdata, self->numCells[0],self->numCells[1],self->numCells[2], self->numCells[0]*self->numCells[1]*self->numCells[2] );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 8, expected key 'POINT_DATA', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==9){
         if (!strcmp(identifier, "SCALARS")){
            if(!self->datasetName)  /** if the user does not specify a dataset, use the first one */
               self->datasetName = StG_Strdup( strtok_r(NULL, self->delim, &self->currentStrtokKey ) );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 9, expected key 'SCALARS', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }
      if(jj==10){
         if (!strcmp(identifier, "LOOKUP_TABLE")){
            if ( strcmp( "default", strtok_r(NULL, self->delim, &self->currentStrtokKey ) ) )
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nError encountered for vtk file '%s'\n"
                                                       "Lookup table must be of the type 'default'.\n ", \
                                                       __func__,self->type,self->name, self->fileMeta->name );
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                                    "At line 10, expected key 'LOOKUP_TABLE', but instead encountered '%s'\n "
                                                    STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileMeta->name, identifier );
      }

   }

   /** vtk files specify the origin for the domain, whereas our internal book keeping startCrd corresponds to the first point's position. adjust accordingly */
   self->startCrd[0] += 0.5*self->voxelSize[0];
   self->startCrd[1] += 0.5*self->voxelSize[1];
   self->startCrd[2] += 0.5*self->voxelSize[2];

   return 1;
}

int _VoxelDataHandler_VTKStructuredPoints_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_VTKStructuredPoints*  self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;
   char* result;
   char* token;
   char* foundDatasets=NULL;
   int done = 0;
   int doneInner = 0;

   /** if open, close, then reopen file */
   if( self->fileData) 
      File_Close( self->fileData );
   self->fileData  = CFile_NewRead( self->filename );

   /** keep reading until we encounter the lookup table line.  the next line will contain the data */
   while (done==0){
      do {
         result = fgets( self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
         if( self->currentString != result ){
            if( !foundDatasets ) 
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s' before any datasets could be read.\n"
                                                    "Please ensure that the file is not corrupted.",__func__,self->type,self->name, self->fileData->name );
            else
               Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnexpectedly reached end of data file '%s'.\n"
                                                       "The requested dataset was '%s', but only these datasets were found:\n\n%s\n"
                                                       "Please ensure that the correct dataset was specifed (via the DatasetName parameter), or that file is not corrupted.",__func__,self->type,self->name, self->fileData->name, self->datasetName, foundDatasets );               
         }
      token = strtok_r(self->currentString, self->delim, &self->currentStrtokKey );
      if(token)
         doneInner = !strcmp( token, "SCALARS" );
      } while ( !doneInner );

      token = strtok_r(NULL, self->delim, &self->currentStrtokKey );
      if( !strcmp( token, self->datasetName ) ){
         done = 1;
         /* if required dataset, proceed to determine other info */
         char* dtype = strtok_r(NULL, self->delim, &self->currentStrtokKey );
         if(              !strcmp(dtype, "char") ){
            self->dataType = Variable_DataType_Char;
         } else if(       !strcmp(dtype, "unsigned_char") ){
            self->dataType = Variable_DataType_Char;
            self->unsignedType = 1;
         } else if(       !strcmp(dtype, "int") ){
            self->dataType = Variable_DataType_Int;
         } else if(       !strcmp(dtype, "float") ){
            self->dataType = Variable_DataType_Float;
         } else if(       !strcmp(dtype, "double") ){
            self->dataType = Variable_DataType_Double;
         } else
            Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nError encountered for vtk file '%s'\n"
                                                    "The data format '%s' is not supported.\n", \
                                                    __func__,self->type,self->name, self->fileData->name, dtype);

      } else {
         if( !foundDatasets)
            Stg_asprintf( &foundDatasets, "" );

         Stg_asprintf( &foundDatasets, "%s%s\n", foundDatasets, token );
      }
   }

   if (foundDatasets) Memory_Free(foundDatasets);


   /** now read one more line, and check that it is as expected */
   fgets( self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );
   token = strtok_r(self->currentString, self->delim, &self->currentStrtokKey );
   if (!strcmp(token, "LOOKUP_TABLE")){
      if ( strcmp( "default", strtok_r(NULL, self->delim, &self->currentStrtokKey ) ) )
         Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nError encountered for vtk file '%s'\n"
                                                 "Lookup table for dataset '%s' must be of the type 'default'.\n ", \
                                                 __func__,self->type,self->name, self->fileData->name, self->datasetName );
   } else
      Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s'\nUnable to read meta data from file '%s'\n"
                                              "Expected key 'LOOKUP_TABLE' for dataset '%s', but instead encountered '%s'\n "
                                              STANDARDVTKFILECOMMENT ,__func__,self->type,self->name, self->fileData->name, self->datasetName, token );

   /* init cursor */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( voxelDataHandler );

   return 1;
}

int _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataASCII( void* voxelDataHandler ){
   VoxelDataHandler_VTKStructuredPoints*  self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;

   _VoxelDataHandler_Abstract_GetASCIIVoxelData( voxelDataHandler );

   switch ( self->dataType )
   {
      case Variable_DataType_Char   :sscanf( self->currentToken,  "%c", (  char*)self->dataBuffer ); break;
      case Variable_DataType_Int    :sscanf( self->currentToken,  "%i", (   int*)self->dataBuffer ); break;
      case Variable_DataType_Float  :sscanf( self->currentToken,  "%f", ( float*)self->dataBuffer ); break;
      case Variable_DataType_Double :sscanf( self->currentToken, "%lf", (double*)self->dataBuffer ); break;
   }

   return 1;
}

int _VoxelDataHandler_VTKStructuredPoints_GetCurrentVoxelDataBinary( void* voxelDataHandler ){
   VoxelDataHandler_VTKStructuredPoints*  self = (VoxelDataHandler_VTKStructuredPoints*)voxelDataHandler;

   _VoxelDataHandler_Abstract_GetBinaryVoxelData( voxelDataHandler );

   if(!self->unsignedType)
      memcpy( self->dataBuffer, self->tempDataBuffer, self->tempDataBufferSize );
   else{
      /** only currently support unsigned char */
      *(char*)self->dataBuffer = (char) *(unsigned char*)self->tempDataBuffer;
   }

   return 1;
}

void _VoxelDataHandler_VTKStructuredPoints_SkipBlank( File* file ){
   int                              testChar;
   Bool                             blankLine = True;
   /** first check if line is empty.  keep skipping until non-empty line is found */
   while( blankLine == True ){
      testChar = fgetc( CFile_Ptr( file ) );
      if ( testChar != '\n' )
         blankLine = False;
   }
   /** we must now rewind a character */
   fseek( CFile_Ptr( file ), -1, SEEK_CUR );
}
