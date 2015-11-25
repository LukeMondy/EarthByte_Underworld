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
#include "VoxelDataHandler_GocadAbstract.h"
#include "VoxelDataHandler_GocadProperties.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_GocadProperties_Type = "VoxelDataHandler_GocadProperties";

VoxelDataHandler_GocadProperties* _VoxelDataHandler_GocadProperties_New( VOXELDATAHANDLER_GOCADPROPERTIES_DEFARGS )
{
   VoxelDataHandler_GocadProperties* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_GocadProperties ) );
   self = (VoxelDataHandler_GocadProperties*)_VoxelDataHandler_GocadAbstract_New(  VOXELDATAHANDLER_GOCADABSTRACT_PASSARGS  );

   self->requiredPropertyName = NULL;
   sscanf( "Float", "%s",  self->property_datatype ); /** default to Float if nothing is specified */
   self->bigEndian          = 1;
   return self;
}


void* _VoxelDataHandler_GocadProperties_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_GocadProperties);
   Type                                                                        type = VoxelDataHandler_GocadProperties_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_GocadAbstract_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_GocadAbstract_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_GocadProperties_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_GocadProperties_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_GocadAbstract_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_GocadAbstract_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_GocadAbstract_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_GocadProperties_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_GocadAbstract_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_GocadProperties_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_GocadProperties_New( VOXELDATAHANDLER_GOCADPROPERTIES_PASSARGS  );
}

void _VoxelDataHandler_GocadProperties_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_GocadProperties*  self = (VoxelDataHandler_GocadProperties*) voxelDataHandler;
   char* requiredPropertyName;
   unsigned requiredPropertyID;

   /** config parent */
   _VoxelDataHandler_GocadAbstract_AssignFromXML( voxelDataHandler, cf, data );
   /** The name of the required property.  If none is provided, we will use the property with ID=1 */
   requiredPropertyName = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"PropertyName", NULL ));
   /** Alternatively the ID of the required property may be provided.  Note that if a name is provided,         */
   /** this will override any provided property ID.  If neither is provided, we set the default propertyID to 1 */
   requiredPropertyID = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"PropertyID", 1 );

   _VoxelDataHandler_GocadProperties_Init( self, requiredPropertyName, requiredPropertyID );

}

void _VoxelDataHandler_GocadProperties_Init( VoxelDataHandler_GocadProperties* voxelDataHandler, char* requiredPropertyName, unsigned requiredPropertyID ){
   VoxelDataHandler_GocadProperties*  self = (VoxelDataHandler_GocadProperties*) voxelDataHandler;

   self->requiredPropertyName = requiredPropertyName;
   self->requiredPropertyID   = requiredPropertyID;
   /** go ahead and load meta data */
   _VoxelDataHandler_GocadProperties_LoadMetaData( voxelDataHandler );

   /** now open binary/heavy data file */
   self->fileData  = CFile_NewRead( self->filenameData );
   Journal_Firewall( self->fileData != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open data file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filenameData );

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
   self->dataBuffer = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "voxel data handler type" );;

}

int _VoxelDataHandler_GocadProperties_LoadMetaData( void* voxelDataHandler ){
   VoxelDataHandler_GocadProperties*  self = (VoxelDataHandler_GocadProperties*)voxelDataHandler;
   unsigned esize;
   char etype[20];
   char currentString[1000];
   char filenametemp[100];
   char prop_format[20];
   char* identifier;

   /** open file */
   if( self->fileMeta ) {
      File_Reopen( self->fileMeta );
   }
   else {
      self->fileMeta  = CFile_NewRead( self->filename );
   }
   Journal_Firewall( self->fileMeta != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open header file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filename );

   /** if a property name is provided, search for this property and finds its ID */
   if(self->requiredPropertyName){
      HashTable* hashTbl = HashTable_New( NULL, NULL, NULL, HASHTABLE_STRING_KEY );
      HashTable_Index* currentHash;
      unsigned* hashKeyData;
      /** lets first find all the property datasets available */
      /** rewind file to top */
      rewind( CFile_Ptr( self->fileMeta ) );
      /** now read in data from file */
      /** tokenize string */
      while( fgets( currentString, 1000, CFile_Ptr( self->fileMeta ) ) ) {
         /** get row identifier */
         identifier = strtok(currentString, "\" " );
         if(!strcmp(identifier, "PROPERTY")){
            unsigned propertyID;
            void* dataPtr = malloc(sizeof(unsigned));
            char* namePtr;
            sscanf( strtok(NULL, "\" "), "%u" , &propertyID );                     /** sscanf to convert string to unsigned int */
            namePtr = strtok(NULL, "\" ");                                         /** copy pointer to next token, which will contain name of property */
            dataPtr = memcpy(dataPtr,(const void*)&propertyID,sizeof(unsigned));   /** create copy of data which will be stored within hash table entry */
            HashTable_InsertEntry( hashTbl, namePtr, strlen(namePtr), dataPtr, sizeof(unsigned));
         }
      }
      /** now search for entry with provided name */
      hashKeyData = HashTable_FindEntry( hashTbl, self->requiredPropertyName, strlen(self->requiredPropertyName), unsigned);
      if( hashKeyData )
         self->requiredPropertyID = *hashKeyData;
      else {
         currentHash = HashTable_First( hashTbl );
         if(currentHash)
            printf("\n\nProvided property name %s not found.  Properties available in voxel file %s are:\n\n", self->requiredPropertyName, self->filename);
         else
            printf("\n\nProvided property name %s not found.  No properties found in voxel file %s.\nProperties are signified by the string PROPERTY.\n", self->requiredPropertyName, self->filename);
         while( currentHash != NULL ){
            printf("Property Id = %3d, Property Name = %s\n", *(unsigned*)currentHash->curr->data, (char*)currentHash->curr->key);
            currentHash = HashTable_Next( currentHash );
         }
         printf("\n");
         Stg_Class_Delete( hashTbl );
         Journal_Firewall( 0,
            self->errorStream,
            "Error in %s for %s '%s'.",
            __func__,
            self->type,
            self->name );
      }
      Stg_Class_Delete( hashTbl );
   }

   /** now get the information for the required property */
   /** rewind file to top */
   rewind( CFile_Ptr( self->fileMeta ) );
   /** now read in data from file */
   /** tokenize string */
   while( fgets( currentString, 1000, CFile_Ptr( self->fileMeta ) ) ) {
      unsigned currentID;
      /** get row identifier */
      identifier = strtok(currentString, " " );
      if(      !strcmp(identifier, "PROPERTY_SUBCLASS")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID){   /** ensure required property is being read */
            strtok(NULL, " "); /** skip token */
            sscanf( strtok(NULL, " "), "%s",  self->property_datatype );
         }
      }else if(!strcmp(identifier, "PROP_NO_DATA_VALUE")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID)   /** ensure required property is being read */
            sscanf( strtok(NULL, " "), "%lg",  &self->nodata_value );
      }else if(!strcmp(identifier, "PROP_ESIZE")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID)   /** ensure required property is being read */
         sscanf( strtok(NULL, " "), "%u",  &esize );
      }else if(!strcmp(identifier, "PROP_ETYPE")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID)   /** ensure required property is being read */
         sscanf( strtok(NULL, " "), "%s",  etype );
      }else if(!strcmp(identifier, "PROP_FORMAT")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID)   /** ensure required property is being read */
         sscanf( strtok(NULL, " "), "%s",  prop_format );
      }else if(!strcmp(identifier, "PROP_FILE")){
         sscanf( strtok(NULL, " "), "%u",  &currentID );
         if(currentID == self->requiredPropertyID)   /** ensure required property is being read */
         sscanf( strtok(NULL, " "), "%s",  &filenametemp );
      }
   }
   sprintf( self->filenameData, "%s", filenametemp);

   if(     !strcmp(self->property_datatype, "Char") ){
      self->dataType = Variable_DataType_Char;
      Journal_Firewall( esize==1,
         self->errorStream,
         "Error in %s for %s '%s' - Datatype char expecting esize 1 but found etype %u.",
         __func__,
         self->type,
         self->name,
         esize);
   } else if(!strcmp(self->property_datatype, "Int"   )){
      self->dataType = Variable_DataType_Int;
      Journal_Firewall( esize==4,
         self->errorStream,
         "Error in %s for %s '%s' - Datatype Int expecting esize 4 but found etype %u.",
         __func__,
         self->type,
         self->name,
         esize);
   } else if(!strcmp(self->property_datatype, "Float" )){
      self->dataType = Variable_DataType_Float;
      Journal_Firewall( esize==4,
         self->errorStream,
         "Error in %s for %s '%s' - Datatype Float expecting esize 4 but found etype %u.",
         __func__,
         self->type,
         self->name,
         esize);
   } else if(!strcmp(self->property_datatype, "Double")){
      self->dataType = Variable_DataType_Double;
      Journal_Firewall( esize==8,
         self->errorStream,
         "Error in %s for %s '%s' - Datatype Double expecting esize 8 but found etype %u.",
         __func__,
         self->type,
         self->name,
         esize);
   } else
      Journal_Firewall( False,
        self->errorStream,
        "Error in %s for %s '%s' - Unknown data type (%s) specified within gocad voxel header file.",
        __func__,
        self->type,
        self->name,
        self->property_datatype);

   Journal_Firewall( !strcmp(etype, "IEEE"),
      self->errorStream,
      "Error in %s for %s '%s' - PROP_ETYPE \"%s\" unknown or unsupported.  PROP_ETYPE must be \"IEEE\".",
      __func__,
      self->type,
      self->name,
      etype);

   Journal_Firewall( !strcmp(prop_format, "RAW"),
      self->errorStream,
      "Error in %s for %s '%s' - PROP_FORMAT \"%s\" unknown or unsupported.  PROP_FORMAT must be \"RAW\".",
      __func__,
      self->type,
      self->name,
      prop_format);

   /** should add some checks here */
   File_Close( self->fileMeta );
   return 1;
}

void _VoxelDataHandler_GocadProperties_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_GocadProperties*  self = (VoxelDataHandler_GocadProperties*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_GocadAbstract_Destroy( self, data );
}

int _VoxelDataHandler_GocadProperties_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_GocadProperties*  self = (VoxelDataHandler_GocadProperties*)voxelDataHandler;
   double doubleData;

   _VoxelDataHandler_Abstract_GetBinaryVoxelData( voxelDataHandler );

   memcpy( self->dataBuffer, self->tempDataBuffer, self->tempDataBufferSize );

   switch ( self->dataType )
   {
      case Variable_DataType_Char   : doubleData = (double) *(signed char*)self->dataBuffer; break;
      case Variable_DataType_Int    : doubleData = (double) *(        int*)self->dataBuffer; break;
      case Variable_DataType_Float  : doubleData = (double) *(      float*)self->dataBuffer; break;
      case Variable_DataType_Double : doubleData = (double) *(     double*)self->dataBuffer; break;
      default                       : Journal_Firewall( 0, self->errorStream, "Error in %s for %s '%s' - Datatype not know or not supported.",__func__,self->type,self->name ); break;
   }
   if(doubleData > self->nodata_value)
      return 1;
   else
      return 0;
}

