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
#include "VoxelDataHandler_GocadAbstract.h"
#include "VoxelDataHandler_GocadMaterials.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_GocadMaterials_Type = "VoxelDataHandler_GocadMaterials";

VoxelDataHandler_GocadMaterials* _VoxelDataHandler_GocadMaterials_New( VOXELDATAHANDLER_GOCADMATERIALS_DEFARGS )
{
   VoxelDataHandler_GocadMaterials* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_GocadMaterials ) );
   self = (VoxelDataHandler_GocadMaterials*)_VoxelDataHandler_GocadAbstract_New(  VOXELDATAHANDLER_GOCADABSTRACT_PASSARGS  );

   self->dataType           = Variable_DataType_Char;
   self->dataBufferSize     = sizeof(signed char);
   self->dataBuffer         = Memory_Alloc_Bytes_Unnamed( self->dataBufferSize , "gocad voxel data handler type" );
   self->numberRegions      = 0;
   self->esize              = 0;
   self->matIndexMapping    = NULL;
   self->hashTbl            = HashTable_New( NULL, NULL, NULL, HASHTABLE_STRING_KEY );
   self->bigEndian          = 0;

   return self;
}


void* _VoxelDataHandler_GocadMaterials_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                                _sizeOfSelf = sizeof(VoxelDataHandler_GocadMaterials);
   Type                                                                        type = VoxelDataHandler_GocadMaterials_Type;
   Stg_Class_DeleteFunction*                                                _delete = _VoxelDataHandler_GocadMaterials_Delete;
   Stg_Class_PrintFunction*                                                  _print = _VoxelDataHandler_GocadAbstract_Print;
   Stg_Class_CopyFunction*                                                    _copy = _Stg_Component_Copy;
   AllocationType                                                nameAllocationType = NON_GLOBAL;
   Stg_Component_DefaultConstructorFunction*                    _defaultConstructor = _VoxelDataHandler_GocadMaterials_DefaultNew;
   Stg_Component_ConstructFunction*                                      _construct = _VoxelDataHandler_GocadMaterials_AssignFromXML;
   Stg_Component_BuildFunction*                                              _build = _VoxelDataHandler_GocadMaterials_Build;
   Stg_Component_InitialiseFunction*                                    _initialise = _VoxelDataHandler_GocadAbstract_Initialise;
   Stg_Component_ExecuteFunction*                                          _execute = _VoxelDataHandler_GocadAbstract_Execute;
   Stg_Component_DestroyFunction*                                          _destroy = _VoxelDataHandler_GocadMaterials_Destroy;
   VoxelDataHandler_GetTotalVoxelCount*                         _getTotalVoxelCount = _VoxelDataHandler_Abstract_GetTotalVoxelCount;
   VoxelDataHandler_GotoVoxelDataTop*                             _gotoVoxelDataTop = _VoxelDataHandler_GocadAbstract_GotoVoxelDataTop;
   VoxelDataHandler_IncrementVoxelIndex*                       _incrementVoxelIndex = _VoxelDataHandler_Abstract_IncrementVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelIndex*                     _getCurrentVoxelIndex = _VoxelDataHandler_Abstract_GetCurrentVoxelIndex;
   VoxelDataHandler_GetCurrentVoxelCoord*                     _getCurrentVoxelCoord = _VoxelDataHandler_Abstract_GetCurrentVoxelCoord;
   VoxelDataHandler_GetCurrentVoxelDataFunc*                   _getCurrentVoxelData = _VoxelDataHandler_GocadMaterials_GetCurrentVoxelData;

   return (void*)_VoxelDataHandler_GocadMaterials_New( VOXELDATAHANDLER_GOCADMATERIALS_PASSARGS  );
}

void _VoxelDataHandler_GocadMaterials_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_GocadMaterials* self = (VoxelDataHandler_GocadMaterials*)voxelDataHandler;

   Stg_Class_Delete( self->hashTbl );
   Memory_Free( self->matIndexMapping );
   /* delete parent class */
   _VoxelDataHandler_GocadAbstract_Delete( self );
}

void _VoxelDataHandler_GocadMaterials_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*) voxelDataHandler;

   /** config parent */
   _VoxelDataHandler_GocadAbstract_AssignFromXML( voxelDataHandler, cf, data );

   _VoxelDataHandler_GocadMaterials_Init( self );

}

void _VoxelDataHandler_GocadMaterials_Init( VoxelDataHandler_GocadMaterials* voxelDataHandler ){
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*) voxelDataHandler;

   /** go ahead and load meta data */
   _VoxelDataHandler_GocadMaterials_LoadMetaData( voxelDataHandler );

   /** now open binary/heavy data file */
   self->fileData  = CFile_NewRead( self->filenameData );
   Journal_Firewall( self->fileData != NULL,
        self->errorStream,
        "Error in %s for %s '%s' - Cannot find or open data file '%s'.",
        __func__,
        self->type,
        self->name,
        self->filenameData );

}

int _VoxelDataHandler_GocadMaterials_LoadMetaData( void* voxelDataHandler ){
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*)voxelDataHandler;
   char currentString[1000];
   char filenametemp[100];
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

   /** tokenize string */
   while( fgets( currentString, 1000, CFile_Ptr( self->fileMeta ) ) ) {
      /** get row identifier */
      identifier = strtok(currentString, " " );
      if(     !strcmp(identifier, "FLAGS_ESIZE"))
         sscanf( strtok(NULL, " "), "%u",  &self->esize );
      else if(!strcmp(identifier, "FLAGS_FILE"))
         sscanf( strtok(NULL, " "), "%s",  &filenametemp );
      else if(!strcmp(identifier, "REGION") || !strcmp(identifier, "MODEL_REGION")){
         unsigned regionID;
         void* dataPtr = malloc(sizeof(unsigned));
         char* namePtr;
         namePtr = strtok(NULL, "\" ");                                       /** copy pointer to region name token */
         sscanf( strtok(NULL, "\" "), "%u" , &regionID );                     /** sscanf to convert string to unsigned int */
         dataPtr = memcpy(dataPtr,(const void*)&regionID,sizeof(unsigned));   /** create copy of data which will be stored within hash table entry */
         HashTable_InsertEntry( self->hashTbl, namePtr, strlen(namePtr), dataPtr, sizeof(unsigned));
      }
   }
   self->numberRegions = self->hashTbl->count;
   sprintf( self->filenameData, "%s", filenametemp);
   Journal_Firewall( self->numberRegions != 0,
        self->errorStream,
        "Error in %s for %s '%s' - No Voxet regions found within file '%s'.\nRegions are signified by the string REGION or MODEL_REGION.",
        __func__,
        self->type,
        self->name,
        self->filename );

   File_Close( self->fileMeta );
   return 1;
}

void _VoxelDataHandler_GocadMaterials_Build( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*)voxelDataHandler;
   HashTable_Index* currentHash;
   Materials_Register*  materials_Register = ((PICelleratorContext*)self->context)->materials_Register;
   Material* material;

   /** build parent */
   _VoxelDataHandler_GocadAbstract_Build( voxelDataHandler, data );

   /** create mapping from position in flags datum to materialIndex */
   self->matIndexMapping = Memory_Alloc_Array( signed char, self->hashTbl->count, "voxeldatahandler_gocadmaterials__matindexmapping");
   currentHash = HashTable_First( self->hashTbl );
   while( currentHash != NULL ){
      material = Materials_Register_GetByName( materials_Register, (char*)currentHash->curr->key );
      if(material)
         self->matIndexMapping[*(unsigned*)currentHash->curr->data-6] = (signed char) material->index;
      else
         self->matIndexMapping[*(unsigned*)currentHash->curr->data-6] = (signed char) NO_MATERIAL_TAG;
      currentHash = HashTable_Next( currentHash );
   }

   /* now that we know esize, build temp buffer */
   self->tempDataBufferSize = self->esize;
   self->tempDataBuffer     = Memory_Alloc_Bytes_Unnamed( self->tempDataBufferSize , "gocad voxel data handler type" );

}

void _VoxelDataHandler_GocadMaterials_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*)voxelDataHandler;
   /** destroy parent */
   _VoxelDataHandler_GocadAbstract_Destroy( self, data );
}

int _VoxelDataHandler_GocadMaterials_GetCurrentVoxelData( void* voxelDataHandler ){
   VoxelDataHandler_GocadMaterials*  self = (VoxelDataHandler_GocadMaterials*)voxelDataHandler;
   unsigned ii,jj;
   unsigned regionCount=0;

   _VoxelDataHandler_Abstract_GetBinaryVoxelData( voxelDataHandler );

   /** gocad stores materials flags as a series of bits representing true or false, though the first 6 bits in the datum
     are reserved for connectivity data.  note that although in gocad it is possible for multiple materials to exist
     in a region, this is not allowed in the underworld definition of materials.  we currently scan datum until we
     find a material corresponding to a UW material, and use that for the definition.  a safer alternative would be to check to ensure only one
     material is defined as true within the datum.   */
   *((signed char*)self->dataBuffer) = -1;  /** default to no material */
   /** scan through required bytes */
   for(ii=0; ii<self->esize; ii++){
     /** scan through each bit of each byte */
     for(jj = 0; jj < 8; jj++){
        /** skip connectivity data in first byte (first 6 bits), and also stop once reached end of relevent data (leaving only padding) */
        if( (ii>0 || jj>5) ){
           if( *((char*)self->tempDataBuffer+ii) & 1<<jj ){
              *((signed char*)self->dataBuffer) = self->matIndexMapping[regionCount];
              if(self->matIndexMapping[regionCount] != -1) regionCount=self->numberRegions+1; /** if a material is found, we are done, so get out.         */
           }                                                                                  /** note that if no UW material corresponding to the gocad   */
           regionCount++;                                                                     /** material is found, we continue to scan current datum     */
        }                                                                                     /** to determine if a correpsonding materials is found later */
        if(regionCount > self->numberRegions) break;                                          /** within the datum                                         */
     }
     if(regionCount > self->numberRegions) break;
  }

  if(*((signed char*)self->dataBuffer) > -1)
    return 1;
  else
    return 0;
}

