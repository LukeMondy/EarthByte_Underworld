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
#include "VoxelDataHandler_Geomodeller.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_Geomodeller_Type = "VoxelDataHandler_Geomodeller";

VoxelDataHandler_Geomodeller* _VoxelDataHandler_Geomodeller_New( VOXELDATAHANDLER_GEOMODELLER_DEFARGS )
{
   VoxelDataHandler_Geomodeller* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_Geomodeller ) );
   self = (VoxelDataHandler_Geomodeller*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->dataType       = Variable_DataType_Char;
   self->dataBufferSize = sizeof(signed char);
   self->dataBuffer     = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "voxel data handler type" );
   self->dataStride     = -1;
   self->dataPos        = 1;

   return self;
}

void _VoxelDataHandler_Geomodeller_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_Geomodeller* self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;

   /* delete parent class */
   _VoxelDataHandler_Abstract_Delete( self );
}

void _VoxelDataHandler_Geomodeller_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_Geomodeller* self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_Geomodeller (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _VoxelDataHandler_Abstract_Print( self, stream );

   /* VoxelDataHandler_Geomodeller */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void* _VoxelDataHandler_Geomodeller_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_Geomodeller);
   Type                                                                        type = VoxelDataHandler_Geomodeller_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_Geomodeller_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_Geomodeller_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_Geomodeller_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_Geomodeller_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_Geomodeller_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_Geomodeller_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_Geomodeller_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_Geomodeller_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_Geomodeller_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_Geomodeller_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_Geomodeller_New( VOXELDATAHANDLER_GEOMODELLER_PASSARGS  );
}

void _VoxelDataHandler_Geomodeller_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_Geomodeller*    self     = (VoxelDataHandler_Geomodeller*) voxelDataHandler;

   /** config parent */
   _VoxelDataHandler_Abstract_AssignFromXML( voxelDataHandler, cf, data );

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
   _VoxelDataHandler_Geomodeller_LoadMetaData( voxelDataHandler );

   /** now that we've loaded required information, setup auxilliary quantities */
   _VoxelDataHandler_Abstract_CompleteSetup( voxelDataHandler );
}

void _VoxelDataHandler_Geomodeller_Build( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Build( voxelDataHandler, data); }

void _VoxelDataHandler_Geomodeller_Initialise( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Initialise( voxelDataHandler, data); }

void _VoxelDataHandler_Geomodeller_Execute( void* voxelDataHandler, void* data )  { _VoxelDataHandler_Abstract_Execute( voxelDataHandler, data); }

void _VoxelDataHandler_Geomodeller_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_Geomodeller*  self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_Abstract_Destroy( self, data );
}

int _VoxelDataHandler_Geomodeller_LoadMetaData( void* voxelDataHandler ){
   VoxelDataHandler_Geomodeller*  self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;
   int jj;
   char* identifier;
   /** rewind file to top */
   rewind( CFile_Ptr( self->fileMeta ) );
   /** skip blank lines (if any) */
   _VoxelDataHandler_Geomodeller_SkipBlank( self->fileMeta );
   /** now read in data from file */
   /** tokenize string */
   for(jj=0; jj<10 ; jj++) {
      /** get string from file */
      fgets( self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileMeta ) );
      /** get row identifier */
      identifier = strtok_r(self->currentString, self->delim, &self->currentStrtokKey );
      if(      !strcmp(identifier, "nx"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u",  &self->numCells[0] );
      else if (!strcmp(identifier, "ny"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u",  &self->numCells[1] );
      else if (!strcmp(identifier, "nz"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%u",  &self->numCells[2] );
      else if (!strcmp(identifier, "x0"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[0] );
      else if (!strcmp(identifier, "y0"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[1] );
      else if (!strcmp(identifier, "z0"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->startCrd[2] );
      else if (!strcmp(identifier, "dx"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[0] );
      else if (!strcmp(identifier, "dy"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[1] );
      else if (!strcmp(identifier, "dz"))           sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ), "%lg", &self->voxelSize[2] );
      else if (!strcmp(identifier, "nodata_value")) sscanf( strtok_r(NULL, self->delim, &self->currentStrtokKey ),  "%s", self->nodata_string );
      else Journal_Firewall( 0,
           self->errorStream,
           "Error in %s for %s '%s' - The string \'%s\' within Geomodeller Voxel file %s is not recognised.\n The first 10 lines of your file must specifiy all of nx,ny,nz,x0,y0,z0,dx,dy,dz,nodata_value and nothing else",
           __func__,
           self->type,
           self->name,
           identifier,
           self->filename);
   }

   return 1;
}

int _VoxelDataHandler_Geomodeller_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_Geomodeller*  self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;
   int jj;

   /** if open, close, then reopen file */
   if( self->fileData) 
      File_Close( self->fileData );
   self->fileData  = CFile_NewRead( self->filename );

   rewind( CFile_Ptr( self->fileData ) );
   /** skip blank lines (if any) */
   _VoxelDataHandler_Geomodeller_SkipBlank( self->fileData );
   /** skip meta data section */
   for(jj=0; jj<10 ; jj++)
      fgets( self->currentString, VOXELMAX_LINE_LENGTH_DEFINE, CFile_Ptr( self->fileData ) );

   /* init cursor */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( voxelDataHandler );

   return 1;
}

int _VoxelDataHandler_Geomodeller_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_Geomodeller*  self = (VoxelDataHandler_Geomodeller*)voxelDataHandler;
   Materials_Register*  materials_Register = ((PICelleratorContext*)self->context)->materials_Register;
   Material* material;
   signed char matIndex;

   _VoxelDataHandler_Abstract_GetASCIIVoxelData( voxelDataHandler );

   material = Materials_Register_GetByName( materials_Register, self->currentToken );
   if(material){
      matIndex = (signed char) material->index;
      /** copy data into buffer space */
      memcpy(self->dataBuffer, (void*)&matIndex, sizeof(signed char));
      return 1;
   } else if( !strcmp(self->currentToken, self->nodata_string) ){  /** set material index to -1 if nodata_string is encountered */
      matIndex = (signed char) NO_MATERIAL_TAG;
      /** copy data into buffer space */
      memcpy(self->dataBuffer, (void*)&matIndex, sizeof(signed char));
      return 0;
   } else {
      Journal_Firewall( (long)material,
           self->errorStream,
           "Error in %s for %s '%s' - No material of the name '%s' found.\nAll materials in voxel file must correspond to a StGermain material.",
           __func__,
           self->type,
           self->name,
           self->currentToken );
      return -1; /* should never get here */
    }
}

void _VoxelDataHandler_Geomodeller_SkipBlank( File* file ){
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

