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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

const Type VoxelDataHandler_GocadAbstract_Type = "VoxelDataHandler_GocadAbstract";

VoxelDataHandler_GocadAbstract* _VoxelDataHandler_GocadAbstract_New( VOXELDATAHANDLER_GOCADABSTRACT_DEFARGS )
{
   VoxelDataHandler_GocadAbstract* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof( VoxelDataHandler_GocadAbstract ) );
   self = (VoxelDataHandler_GocadAbstract*)_VoxelDataHandler_Abstract_New(  VOXELDATAHANDLER_ABSTRACT_PASSARGS  );

   self->cursorData   = 0;
   self->currentData  = NULL;
   self->zPositiveUp  = True;
   /** organise default coords */
   self->coord_axis[0][0] = 1;
   self->coord_axis[0][1] = 0;
   self->coord_axis[0][2] = 0;
   self->coord_axis[1][0] = 0;
   self->coord_axis[1][1] = 1;
   self->coord_axis[1][2] = 0;
   self->coord_axis[2][0] = 0;
   self->coord_axis[2][1] = 0;
   self->coord_axis[2][2] = 1;

   /** init numcells so we can tell if they've ever been set */
   self->numCells[0] = (unsigned) -1;
   self->numCells[1] = (unsigned) -1;
   self->numCells[2] = (unsigned) -1;
   /** init cell size so we can tell if they've ever been set */
   self->voxelSize[0] = -1;
   self->voxelSize[1] = -1;
   self->voxelSize[2] = -1;

   return self;
}

void _VoxelDataHandler_GocadAbstract_Delete( void* voxelDataHandler ) {
   VoxelDataHandler_GocadAbstract* self = (VoxelDataHandler_GocadAbstract*)voxelDataHandler;

   /* delete parent class */
   _VoxelDataHandler_Abstract_Delete( self );
}

void _VoxelDataHandler_GocadAbstract_Print( void* voxelDataHandler, Stream* stream ) {
   VoxelDataHandler_GocadAbstract* self = (VoxelDataHandler_GocadAbstract*)voxelDataHandler;

   /* General info */
   Journal_Printf( stream, "VoxelDataHandler_GocadAbstract (ptr): %p:\n", self );
   Stream_Indent( stream );

   /* Parent class info */
   _VoxelDataHandler_Abstract_Print( self, stream );

   /* VoxelDataHandler_GocadAbstract */
   Journal_Printf( stream, "filename: %s\n", self->filename );

   Stream_UnIndent( stream );
}

void _VoxelDataHandler_GocadAbstract_AssignFromXML( void* voxelDataHandler, Stg_ComponentFactory *cf, void* data ) {
   VoxelDataHandler_GocadAbstract*  self = (VoxelDataHandler_GocadAbstract*) voxelDataHandler;

   /** config parent */
   _VoxelDataHandler_Abstract_AssignFromXML( voxelDataHandler, cf, data );

   _VoxelDataHandler_GocadAbstract_Init( self );

}

void _VoxelDataHandler_GocadAbstract_Init( VoxelDataHandler_GocadAbstract* voxelDataHandler ){
   VoxelDataHandler_GocadAbstract*  self = (VoxelDataHandler_GocadAbstract*) voxelDataHandler;
   unsigned ii, jj;

   /** go ahead and load meta data */
   _VoxelDataHandler_GocadAbstract_LoadMetaData( voxelDataHandler );
   for(ii=0; ii<3; ii++)
      for(jj=0; jj<3; jj++) Journal_Firewall( !((abs(self->coord_axis[ii][jj]) && (ii!=jj)) ) ,
                                self->errorStream,
                                "Error in %s for %s '%s' - Currently general voxel coordinates are not supported.  UVW axis must respectively be parallel to XYZ coordinates.",
                                __func__,
                                self->type );

   /** determine voxelsize or numcells (whichever is missing ) */
   if( (self->voxelSize[0] <= 0) || (self->voxelSize[1] <= 0) || (self->voxelSize[2] <= 0) )
      for(ii=0; ii<3; ii++)  self->voxelSize[ii] = ( self->axis_max[ii] - self->axis_min[ii] )*self->coord_axis[ii][ii] / ( self->numCells[ii] - 1 );
   else if ( (self->numCells[0] == (unsigned) -1) || (self->numCells[1] == (unsigned) -1) || (self->numCells[2] == (unsigned) -1) )
      for(ii=0; ii<3; ii++)  self->numCells[ii] = ( self->axis_max[ii] - self->axis_min[ii] )*self->coord_axis[ii][ii]  / self->voxelSize[ii]  +  1;
   else
      Journal_Firewall( False,
        self->errorStream,
        "Error in %s for %s '%s' - Either AXIS_N or AXIS_D must be specified within gocad voxel header file.",
        __func__,
        self->type,
        self->name );

   /** note that the following is also performed in _VoxelDataHandler_Abstract_CalcMinMaxCoords, but do it here as well just to be clear */
   for(ii=0; ii<3; ii++){
      /** ensure voxelSize is positive, and juggle things if it isn't */
      if(self->voxelSize[ii]<0){
         self->scale[ii]     = -self->scale[ii]; /** reverse the scale so that we count in opposite direction */
         self->voxelSize[ii] = -self->voxelSize[ii]; /** set voxelSize to be positive */
      }
   }

   /** if zPositiveUp is false (ie, z axis is positive downwards), negative the scale to reflect this */
   if(self->zPositiveUp == False) self->scale[2] = -self->scale[2];

   /** organise start crds.  note that gocad voxel data is specified on the node (ie, cell centres), so no need to adjust */
   for(ii=0; ii<3; ii++) self->startCrd[ii] +=  self->axis_min[ii]*self->coord_axis[ii][ii];

   /** complete setup! */
   _VoxelDataHandler_Abstract_CompleteSetup( voxelDataHandler );

}

int _VoxelDataHandler_GocadAbstract_LoadMetaData( void* voxelDataHandler ){
   VoxelDataHandler_GocadAbstract*  self = (VoxelDataHandler_GocadAbstract*)voxelDataHandler;
   char currentString[1000];
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

   /** rewind file to top */
   rewind( CFile_Ptr( self->fileMeta ) );
   /** now read in data from file */
   /** tokenize string */
   while( fgets( currentString, 1000, CFile_Ptr( self->fileMeta ) ) ) {
      /** get row identifier */
      identifier = strtok(currentString, " " );
      if(      !strcmp(identifier, "ZPOSITIVE")){
         if(!strcmp(strtok(NULL, " "), "Depth"))
            self->zPositiveUp = False;
         else
            self->zPositiveUp = True;
      }else if(!strcmp(identifier, "AXIS_O")){
         sscanf( strtok(NULL, " "), "%lg",  &self->startCrd[0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->startCrd[1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->startCrd[2] );
      }else if(!strcmp(identifier, "AXIS_U")){
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[0][0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[0][1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[0][2] );
      }else if(!strcmp(identifier, "AXIS_V")){
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[1][0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[1][1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[1][2] );
      }else if(!strcmp(identifier, "AXIS_W")){
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[2][0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[2][1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->coord_axis[2][2] );
      }else if(!strcmp(identifier, "AXIS_N")){
         sscanf( strtok(NULL, " "), "%u",  &self->numCells[0] );
         sscanf( strtok(NULL, " "), "%u",  &self->numCells[1] );
         sscanf( strtok(NULL, " "), "%u",  &self->numCells[2] );
      }else if(!strcmp(identifier, "AXIS_D")){
         sscanf( strtok(NULL, " "), "%lg",  &self->voxelSize[0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->voxelSize[1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->voxelSize[2] );
      }else if(!strcmp(identifier, "AXIS_MIN")){
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_min[0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_min[1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_min[2] );
      }else if(!strcmp(identifier, "AXIS_MAX")){
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_max[0] );
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_max[1] );
         sscanf( strtok(NULL, " "), "%lg",  &self->axis_max[2] );
      }
   }

   File_Close( self->fileMeta );
   return 1;
}

void _VoxelDataHandler_GocadAbstract_Build( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Build( voxelDataHandler, data); }

void _VoxelDataHandler_GocadAbstract_Initialise( void* voxelDataHandler, void* data ) { _VoxelDataHandler_Abstract_Initialise( voxelDataHandler, data); }

void _VoxelDataHandler_GocadAbstract_Execute( void* voxelDataHandler, void* data )  { _VoxelDataHandler_Abstract_Execute( voxelDataHandler, data); }

void _VoxelDataHandler_GocadAbstract_Destroy( void* voxelDataHandler, void* data ) {
   VoxelDataHandler_GocadAbstract*  self = (VoxelDataHandler_GocadAbstract*)voxelDataHandler;

   /** close file */
   if( self->fileData ) {
      Stg_Class_Delete( self->fileData );
      self->fileData = NULL;
   }

   /** destroy parent */
   _VoxelDataHandler_Abstract_Destroy( voxelDataHandler, data );
}

int _VoxelDataHandler_GocadAbstract_GotoVoxelDataTop( void* voxelDataHandler ){
   VoxelDataHandler_GocadAbstract*  self = (VoxelDataHandler_GocadAbstract*)voxelDataHandler;

   /** if open, close, then reopen file */
   if( self->fileData ) 
      File_Close( self->fileData );

   self->fileData  = CFile_NewRead( self->filenameData );

   rewind( CFile_Ptr( self->fileData ) );

   /* init cursor */
   _VoxelDataHandler_Abstract_GotoVoxelDataTop( voxelDataHandler );
   return 1;
}

