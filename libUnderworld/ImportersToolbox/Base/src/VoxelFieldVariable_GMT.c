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
#include "VoxelFieldVariable.h"
#include "VoxelFieldVariable_GMT.h"

#include <assert.h>
#include <string.h>

const Type VoxelFieldVariable_GMT_Type = "VoxelFieldVariable_GMT";

void* _VoxelFieldVariable_GMT_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                         _sizeOfSelf = sizeof(VoxelFieldVariable_GMT);
   Type                                                                 type = VoxelFieldVariable_GMT_Type;
   Stg_Class_DeleteFunction*                                         _delete = _FieldVariable_Delete;
   Stg_Class_PrintFunction*                                           _print = _FieldVariable_Print;
   Stg_Class_CopyFunction*                                             _copy = _FieldVariable_Copy;
   Stg_Component_DefaultConstructorFunction*             _defaultConstructor = _VoxelFieldVariable_DefaultNew;
   Stg_Component_ConstructFunction*                               _construct = _VoxelFieldVariable_GMT_AssignFromXML;
   Stg_Component_BuildFunction*                                       _build = _VoxelFieldVariable_Build;
   Stg_Component_InitialiseFunction*                             _initialise = _VoxelFieldVariable_Initialise;
   Stg_Component_ExecuteFunction*                                   _execute = _VoxelFieldVariable_Execute;
   Stg_Component_DestroyFunction*                                   _destroy = _VoxelFieldVariable_Destroy;
   AllocationType                                         nameAllocationType = NON_GLOBAL;
   FieldVariable_InterpolateValueAtFunction*             _interpolateValueAt = _VoxelFieldVariable_GMT_InterpolateValueAt;/** wrap this one */
   FieldVariable_GetValueFunction*               _getMinGlobalFieldMagnitude = _VoxelFieldVariable_GetMinGlobalFieldMagnitude;
   FieldVariable_GetValueFunction*               _getMaxGlobalFieldMagnitude = _VoxelFieldVariable_GetMaxGlobalFieldMagnitude;
   FieldVariable_CacheValuesFunction*       _cacheMinMaxGlobalFieldMagnitude = _VoxelFieldVariable_CacheMinMaxGlobalFieldMagnitude;
   FieldVariable_GetCoordFunction*                  _getMinAndMaxLocalCoords = _VoxelFieldVariable_GetMinAndMaxLocalCoords;
   FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords = _VoxelFieldVariable_GetMinAndMaxGlobalCoords;
   VoxelFieldVariable_SetupDataArray*                    _setupDataArray = _VoxelFieldVariable_SetupDataArray;
   VoxelFieldVariable_TestDataArray*                      _testDataArray = _VoxelFieldVariable_TestDataArray;
   VoxelFieldVariable_SetArrayValue*                      _setArrayValue = _VoxelFieldVariable_SetArrayValue;
   VoxelFieldVariable_GetArrayIndexIJK*                _getArrayIndexIJK = _VoxelFieldVariable_GetArrayIndexIJK;
   
   return _VoxelFieldVariable_GMT_New(  VOXELFIELDVARIABLE_GMT_PASSARGS  );
}

VoxelFieldVariable_GMT* _VoxelFieldVariable_GMT_New(  VOXELFIELDVARIABLE_GMT_DEFARGS  ) {
   VoxelFieldVariable_GMT* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(VoxelFieldVariable_GMT) );
   self = (VoxelFieldVariable_GMT*)_VoxelFieldVariable_New(  VOXELFIELDVARIABLE_PASSARGS  );

   self->_setupDataArray   = _setupDataArray;
   self->_testDataArray    = _testDataArray;
   self->_setArrayValue    = _setArrayValue;
   self->_getArrayIndexIJK = _getArrayIndexIJK;
   
   self->localMax = -HUGE_VAL;
   self->localMin =  HUGE_VAL;

   self->decompositionMesh = NULL;

   return self;
}

void _VoxelFieldVariable_GMT_AssignFromXML( void* _VoxelFieldVariable_GMT, Stg_ComponentFactory* cf, void* data ) {
   VoxelFieldVariable_GMT*         self = (VoxelFieldVariable_GMT*)_VoxelFieldVariable_GMT;
   double maxR, maxY;

   _VoxelFieldVariable_AssignFromXML( self, cf, data );
   self->cartesianMode = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"CartesianMode", False );
   self->mapTo360      = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"MapTo360", False );
   self->plateDepth    = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"PlateDepth", 1.0  );

   maxR  =  Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"maxX", 1.0  );
   maxY  =  Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"maxY", 1.0  );

   if(!self->cartesianMode){
       self->innerRadius   = maxR - self->plateDepth;
   }
   else{
       self->innerRadius   = maxY - self->plateDepth;
   }
}

InterpolationResult _VoxelFieldVariable_GMT_InterpolateValueAt( void* voxelFieldVariable, double* globalCoord, double* value ) {
   VoxelFieldVariable_GMT* self = (VoxelFieldVariable_GMT*)voxelFieldVariable;
   InterpolationResult  retValue = OUTSIDE_GLOBAL;
   Coord                coord;
   double               lat,lon,radius;
   double               rads = 180.0/M_PI; /** convert from radians to degrees */
   Bool                 cartesian = self->cartesianMode;

   coord[0] = globalCoord[0];
   coord[1] = globalCoord[1];
   if(self->dim == 2){
       coord[2] = self->zSliceCoord;
       radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
   }
   else{
      coord[2] = globalCoord[2];
      if(!cartesian){ radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2]); }
   }

   coord[2] = 0.0; /* I think for this we always want z = 0 */

   if(!cartesian){
       /* we want to convert coordinates to lat/lon coordinates */
       coord[1] = lat = asin( globalCoord[1]/radius )*rads;  /* Global y coordinate */
       coord[0] = lon = -atan2( globalCoord[2], globalCoord[0] )*rads;/* Global z and x coordinates. Should check for x,y,z = 0: should never happen. */
       if(self->mapTo360){
           coord[0] += 180;
       }
       //printf("coord[0]=%lf coord[1]=%lf\n",coord[0],coord[1]);
   }
   else{
       /** just for moment. to test. let's use x-y as lat lon itself */
       coord[1]=globalCoord[2]; /* latitude */
       coord[0]=globalCoord[0]; /* longitude */
   }
   if(!cartesian){
       if(radius > self->innerRadius){
           retValue = _VoxelFieldVariable_InterpolateValueAt( self, coord, value );
       }
   }
   else{
       if(globalCoord[1] > self->innerRadius){
           retValue = _VoxelFieldVariable_InterpolateValueAt( self, coord, value );
       }       
   }
   return retValue;
}
