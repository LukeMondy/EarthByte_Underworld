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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_ASCII_Type = "VoxelDataHandler_ASCII";

VoxelDataHandler_ASCII* _VoxelDataHandler_ASCII_New( VOXELDATAHANDLER_ASCII_DEFARGS )
{
   VoxelDataHandler_ASCII* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_ASCII ) );
   self = (VoxelDataHandler_ASCII*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->dataType         = Variable_DataType_Double;

   return self;
}

void _VoxelDataHandler_ASCII_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_ASCII* self = (VoxelDataHandler_ASCII*)voxelDataHandler;

   /* delete parent class */
   _VoxelDataHandler_Abstract_Delete( self );
}

void _VoxelDataHandler_ASCII_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_ASCII* self = (VoxelDataHandler_ASCII*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_ASCII (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _VoxelDataHandler_Abstract_Print( self, stream );

   /* VoxelDataHandler_ASCII */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void* _VoxelDataHandler_ASCII_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_ASCII);
   Type                                                                        type = VoxelDataHandler_ASCII_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_ASCII_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_ASCII_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_ASCII_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_ASCII_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_ASCII_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_ASCII_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_ASCII_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_ASCII_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_ASCII_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_ASCII_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_ASCII_New( VOXELDATAHANDLER_ASCII_PASSARGS  );
}

void _VoxelDataHandler_ASCII_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_ASCII*    self     = (VoxelDataHandler_ASCII*) voxelDataHandler;
   char* dataType;

   /** config parent */
   _VoxelDataHandler_Abstract_AssignFromXML( voxelDataHandler, cf, data );
   /** deliminators seperating ascii row data */
   /** broken due to dictionary stripping or escaping characters from xml input.. disable for now */
   //sscanf( Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"fileDeliminator", " ,\t" ), "%[^\n]", self->delim );

   /** see VoxelDataHandler_Abstract header file for info, if not self evident */
   self->numCells[0]  = Stg_ComponentFactory_GetInt(    cf, self->name, (Dictionary_Entry_Key)"NumCellsI", 0 );
   self->numCells[1]  = Stg_ComponentFactory_GetInt(    cf, self->name, (Dictionary_Entry_Key)"NumCellsJ", 0 );
   self->numCells[2]  = Stg_ComponentFactory_GetInt(    cf, self->name, (Dictionary_Entry_Key)"NumCellsK", 0 );
   self->startCrd[0]  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StartCoordI", 0 );
   self->startCrd[1]  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StartCoordJ", 0 );
   self->startCrd[2]  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"StartCoordK", 0 );
   self->voxelSize[0] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"CellSizeI", 0 );
   self->voxelSize[1] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"CellSizeJ", 0 );
   self->voxelSize[2] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"CellSizeK", 0 );
   self->noDataValue  = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"NoDataValue", 0 );

   dataType     = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"DataType", "double" );
   if(     !strcmp(dataType, "char"  )) self->dataType = Variable_DataType_Char;
   else if(!strcmp(dataType, "int"   )) self->dataType = Variable_DataType_Int;
   else if(!strcmp(dataType, "float" )) self->dataType = Variable_DataType_Float;
   else if(!strcmp(dataType, "double")) self->dataType = Variable_DataType_Double;
   else
      Journal_Firewall( NULL, self->errorStream, "Error in %s for %s '%s' - Provided DataType %s not know or not supported.",__func__,self->type,self->name, dataType );
   /** stride for data reading.. so if "PosI1 PosJ1 PosK1 Data1 PosI2 PosJ2 PosK2 Data2...." stride is 4 */
   /** negative values indicate a stride with the allowance of newline as required */
   self->dataStride  = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"DataStride", -1 );
   if( self->dataStride == 0 )
      Journal_Firewall( NULL, self->errorStream, "Error in %s for %s '%s' - A DataStride of zero is not valid. Stride must be >0 or <0",__func__,self->type,self->name);
   /** Position of voxel data in datum. starts counting at 1... so 1,2,3,4..  */
   self->dataPos     = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"DataPos", 1 );
   if( (abs(self->dataStride)<self->dataPos) )
      Journal_Firewall( NULL, self->errorStream, "Error in %s for %s '%s' - Provided DataStride (%i) cannot be less than DataPosition (%i).",__func__,self->type, self->name, self->dataStride, self->dataPos );

   /** now go ahead and complete setup */
   _VoxelDataHandler_Abstract_CompleteSetup( voxelDataHandler );
}

void _VoxelDataHandler_ASCII_Build( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_ASCII*    self     = (VoxelDataHandler_ASCII*) voxelDataHandler;

   _VoxelDataHandler_Abstract_Build( voxelDataHandler, data);

   /** setup voxel data ready to read, incase gototop routine isn't called */
   _VoxelDataHandler_ASCII_GotoVoxelDataTop( self );

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

void _VoxelDataHandler_ASCII_Initialise( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Initialise( voxelDataHandler, data); }

void _VoxelDataHandler_ASCII_Execute( void* voxelDataHandler, void* data )  { _VoxelDataHandler_Abstract_Execute( voxelDataHandler, data); }

void _VoxelDataHandler_ASCII_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_ASCII*  self = (VoxelDataHandler_ASCII*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_Abstract_Destroy( voxelDataHandler, data );
}

int _VoxelDataHandler_ASCII_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_ASCII*  self = (VoxelDataHandler_ASCII*)voxelDataHandler;

   /** if open, close, then reopen file */
   if(self->fileData)
      File_Close( self->fileData );
   self->fileData  = CFile_NewRead( self->filename );
   Journal_Firewall( self->fileData != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filename );

   /** rewind file to top */
   rewind( CFile_Ptr( self->fileData ) );
   /** skip blank lines (if any) */
   _VoxelDataHandler_ASCII_SkipBlank( self->fileData );

   /* now re-init cursor */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( voxelDataHandler );

   return 1;
}

int _VoxelDataHandler_ASCII_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_ASCII*  self = (VoxelDataHandler_ASCII*)voxelDataHandler;
   double doubleData;

   _VoxelDataHandler_Abstract_GetASCIIVoxelData( voxelDataHandler );

   switch ( self->dataType )
   {
      case Variable_DataType_Char   :sscanf( self->currentToken,  "%c", (  char*)self->dataBuffer ); break;
      case Variable_DataType_Int    :sscanf( self->currentToken,  "%i", (   int*)self->dataBuffer ); break;
      case Variable_DataType_Float  :sscanf( self->currentToken,  "%f", ( float*)self->dataBuffer ); break;
      case Variable_DataType_Double :sscanf( self->currentToken, "%lf", (double*)self->dataBuffer ); break;
   }

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

void _VoxelDataHandler_ASCII_SkipBlank( File* file ){
   int   testChar;
   Bool  blankLine = True;
   Bool  headerLine = True;
   char  testStr[100];
   fpos_t strPos;
   int   ii;

   /** first check if line is empty.  keep skipping until non-empty line is found */
   while( blankLine == True ){
      testChar = fgetc(CFile_Ptr( file ));
      if ( testChar != '\n' )
         blankLine = False;
   }
   /** we must now rewind a character */
   fseek( CFile_Ptr( file ), -1, SEEK_CUR );

   /** now skip potential header lines */
   while( headerLine == True ){
      fgetpos(CFile_Ptr( file ), &strPos);
      fgets(&testStr, 100, CFile_Ptr( file ));
      for(ii=0; ii<100; ii++){
         if(      testStr[ii] == ' ') continue;
         else if( testStr[ii] == '#' || testStr[ii] == '>' ) break;
         else{
            headerLine = False;
            fsetpos(CFile_Ptr( file ), &strPos);
            break;
         }
      }
   }
}

