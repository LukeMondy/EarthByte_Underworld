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
#include "VoxelHeightFieldVariable.h"

#include <assert.h>
#include <string.h>

const Type VoxelHeightFieldVariable_Type = "VoxelHeightFieldVariable";

unsigned VHFVmapAXIS[2] = { 0, 2 };

void* _VoxelHeightFieldVariable_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                         _sizeOfSelf = sizeof(VoxelHeightFieldVariable);
   Type                                                                 type = VoxelHeightFieldVariable_Type;
   Stg_Class_DeleteFunction*                                         _delete = _FieldVariable_Delete;
   Stg_Class_PrintFunction*                                           _print = _FieldVariable_Print;
   Stg_Class_CopyFunction*                                             _copy = _FieldVariable_Copy;
   Stg_Component_DefaultConstructorFunction*             _defaultConstructor = _VoxelHeightFieldVariable_DefaultNew;
   Stg_Component_ConstructFunction*                               _construct = _VoxelHeightFieldVariable_AssignFromXML;
   Stg_Component_BuildFunction*                                       _build = _VoxelHeightFieldVariable_Build;
   Stg_Component_InitialiseFunction*                             _initialise = _VoxelHeightFieldVariable_Initialise;
   Stg_Component_ExecuteFunction*                                   _execute = _VoxelHeightFieldVariable_Execute;
   Stg_Component_DestroyFunction*                                   _destroy = _VoxelHeightFieldVariable_Destroy;
   AllocationType                                         nameAllocationType = NON_GLOBAL;
   FieldVariable_InterpolateValueAtFunction*             _interpolateValueAt = _VoxelFieldVariable_InterpolateValueAt;
   FieldVariable_GetValueFunction*               _getMinGlobalFieldMagnitude = _VoxelFieldVariable_GetMinGlobalFieldMagnitude;
   FieldVariable_GetValueFunction*               _getMaxGlobalFieldMagnitude = _VoxelFieldVariable_GetMaxGlobalFieldMagnitude;
   FieldVariable_CacheValuesFunction*       _cacheMinMaxGlobalFieldMagnitude = _VoxelFieldVariable_CacheMinMaxGlobalFieldMagnitude;
   FieldVariable_GetCoordFunction*                  _getMinAndMaxLocalCoords = _VoxelFieldVariable_GetMinAndMaxLocalCoords;
   FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords = _VoxelFieldVariable_GetMinAndMaxGlobalCoords;
   VoxelFieldVariable_SetupDataArray*                        _setupDataArray = _VoxelHeightFieldVariable_SetupDataArray;
   VoxelFieldVariable_TestDataArray*                          _testDataArray = _VoxelHeightFieldVariable_TestDataArray;
   VoxelFieldVariable_SetArrayValue*                          _setArrayValue = _VoxelHeightFieldVariable_SetArrayValue;
   VoxelFieldVariable_GetArrayIndexIJK*                    _getArrayIndexIJK = _VoxelHeightFieldVariable_GetArrayIndexIJK;

   return _VoxelHeightFieldVariable_New(  VOXELHEIGHTFIELDVARIABLE_PASSARGS  );
}

VoxelHeightFieldVariable* _VoxelHeightFieldVariable_New(  VOXELHEIGHTFIELDVARIABLE_DEFARGS  ) {
   VoxelHeightFieldVariable* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(VoxelHeightFieldVariable) );
   self = (VoxelHeightFieldVariable*)_VoxelFieldVariable_New(  VOXELFIELDVARIABLE_PASSARGS  );

   self->voxelDataTestArray = NULL;

   return self;
}


void _VoxelHeightFieldVariable_AssignFromXML( void* _VoxelHeightFieldVariable, Stg_ComponentFactory* cf, void* data ) {
   VoxelHeightFieldVariable*         self = (VoxelHeightFieldVariable*)_VoxelHeightFieldVariable;
   Name                              perpendicularAxis;

   _VoxelFieldVariable_AssignFromXML( self, cf, data );
   /** determine axis perpendicular to the heightfield plane */
   perpendicularAxis = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"perpendicularAxis", "Y" );

   Journal_Firewall(
         perpendicularAxis != NULL && strlen( perpendicularAxis  ) >= 1,
         self->errorStream,
         "\n\Error in func %s for %s '%s' - Bad 'perpendicularAxis' provided.\n\n",
         __func__, self->type, self->name );

   switch ( perpendicularAxis[0] ) {
      case 'X': case 'x':
         self->perpAxis = I_AXIS; break;
      case 'Y': case 'y':
         self->perpAxis = J_AXIS; break;
      case 'Z': case 'z':
         self->perpAxis = K_AXIS; break;
      default:
         Journal_Firewall(
            False,
            self->errorStream,
            "Error in func %s for %s '%s' - Bad 'perpendicularAxis' in dictionary.\n",
            __func__, self->type, self->name );
   }

   Journal_Firewall(
         self->perpAxis == J_AXIS,
         self->errorStream,
         "\n\Error in func %s for %s '%s' - Bad 'perpendicularAxis' provided.\nCurrently only using y-axis as your perpendicular axis is supported\n\n",
         __func__, self->type, self->name );

   self->minimumHeight = (float)Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"MinimumHeight", -HUGE_VAL );
   self->airValueMin   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"AirValueMin", +HUGE_VAL );
   self->airValueMax   = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"AirValueMax", -HUGE_VAL );
   self->topCellOffset = (float)Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"TopCellOffset", -0.5 );
   self->userMinValue  = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"MinimumHeight", False );

   Journal_Firewall( !(self->voxelDataHandler->dataType == Variable_DataType_Double) || !(self->voxelDataHandler->dataType == Variable_DataType_Float), self->errorStream,
                      "\n\Error in func %s for %s '%s' - voxel datatype float/double not supported by this class.\n\n", __func__, self->type, self->name );
   /** as this class stores heigh coordinates, we always uses dataype float to store data */
   self->dataType = Variable_DataType_Float;
   self->dim = 2;  /** dim must be set to two */
}


void _VoxelHeightFieldVariable_SetupDataArray( void* voxelHeightFieldVariable ){
   VoxelHeightFieldVariable* self = (VoxelHeightFieldVariable*)voxelHeightFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   unsigned VHFVmapAXIS[2] = { 0, 2 };
   int ii, jj;

   if( self->decompositionMesh ){
      Journal_Firewall( NULL, self->errorStream, "\n\nError in %s for %s '%s' - While everything is in place, using a decomposition mesh has not been tested! Use with caution.\n\n", __func__, self->type, self->name );
      double displacement[3];
      double meshCrdMin[3], meshCrdMax[3];

      /** first ensure all information is updated if a deformation has occcured */
      Mesh_DeformationUpdate( self->decompositionMesh );
      Mesh_GetLocalCoordRange( self->decompositionMesh, &meshCrdMin, &meshCrdMax );
      for(ii=0; ii<self->dim; ii++){
         /** determine range of voxels required over this processors domain */
         self->arraySize[ii] = (int) ( meshCrdMax[VHFVmapAXIS[ii]] - meshCrdMin[VHFVmapAXIS[ii]] ) / ( abs(self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]])*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[ii]]] ) + 1;
         /** get displacement of mesh from start of voxel data depending on which way we are counting */
         if(self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]] > 0){
            displacement[ii] = meshCrdMin[VHFVmapAXIS[ii]] - self->voxelDataHandler->startCrd[switchAxis[VHFVmapAXIS[ii]]] - 0.5*self->voxelDataHandler->scale[switchAxis[switchAxis[VHFVmapAXIS[ii]]]]*self->voxelDataHandler->voxelSize[switchAxis[switchAxis[VHFVmapAXIS[ii]]]];
            /** determine local support coords */
            self->localCrdMin[ii] = self->voxelDataHandler->minDomainCoordsXYZ[VHFVmapAXIS[ii]] + (int)displacement[ii];
            self->localCrdMax[ii] = self->localCrdMin[ii] + self->arraySize[ii]*abs(self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]])*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[ii]]] ;
            self->localStartCrd[ii] = self->localCrdMin[ii]; 
         } else {
            displacement[ii] = self->voxelDataHandler->startCrd[VHFVmapAXIS[ii]] + 0.5*self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]]*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[ii]]] - meshCrdMax[VHFVmapAXIS[ii]];
            /** determine local support coords */
            self->localCrdMin[ii]   = self->voxelDataHandler->minDomainCoordsXYZ[VHFVmapAXIS[ii]] - (int)displacement[ii];
            self->localCrdMax[ii]   = self->localCrdMax[ii] - self->arraySize[ii]*abs(self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]])*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[ii]]] ;
            self->localStartCrd[ii] = self->localCrdMax[ii]; 
         }
         Journal_Firewall( displacement[ii] > 0, self->errorStream, "\n\nError in %s for %s '%s' - Voxel swarm is defined over a smaller region than local domain, which is not valid.\n\n", __func__, self->type, self->name );
      }
   } else {
      for(ii=0; ii<self->dim; ii++){
         /** local array size equal to the voxel size */
         self->arraySize[ii] = self->voxelDataHandler->numCells[switchAxis[VHFVmapAXIS[ii]]];
         /** local support will be set to the support of the voxel dataset */
         self->localCrdMin[ii] = self->voxelDataHandler->minDomainCoordsXYZ[VHFVmapAXIS[ii]];
         self->localCrdMax[ii] = self->voxelDataHandler->maxDomainCoordsXYZ[VHFVmapAXIS[ii]];
         if(self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]] > 0){
            self->localStartCrd[ii] = self->localCrdMin[ii]; 
         } else {
            self->localStartCrd[ii] = self->localCrdMax[ii];
         }
      }
   }
   /** allocate memory */
   self->voxelDataArray     = Memory_Alloc_3DArray( float, self->arraySize[0], self->arraySize[1], 1, (Name)"voxelheightfield array" );
   self->voxelDataTestArray = Memory_Alloc_3DArray( Bool , self->arraySize[0], self->arraySize[1], 1, (Name)"voxelheightfield test array" );
   /** init memory to  */
   for(ii=0; ii<self->arraySize[0]; ii++){
      for(jj=0; jj<self->arraySize[1]; jj++){
            ((float***)self->voxelDataArray)[ii][jj][0] = self->minimumHeight;  /** init field to provided or default min height */
         ((Bool***)self->voxelDataTestArray)[ii][jj][0] = self->userMinValue;   /** is user min height provided, then no need to test as all points are implicitly set */
      }
   }

}

void _VoxelHeightFieldVariable_SetArrayValue( void* voxelHeightFieldVariable, double* coord ){
   VoxelHeightFieldVariable* self = (VoxelHeightFieldVariable*)voxelHeightFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   int indexIJK[3] = { 0, 0, 0 };
   Coord mappedCoord = { coord[VHFVmapAXIS[0]], coord[VHFVmapAXIS[1]], 0 };
   
   /** if coord is in local range, set value.. note that we need to map the axis, since provided with XYZ coords, but require XZ for VHFV */
   if( self->_getArrayIndexIJK( voxelHeightFieldVariable, &mappedCoord, &indexIJK) ){
      double voxelData;
      Bool   materialExists;
      float* currentHeightValPtr = &((float***)self->voxelDataArray    )[indexIJK[0]][indexIJK[1]][0];
      Bool*    currentTestValPtr = &(( Bool***)self->voxelDataTestArray)[indexIJK[0]][indexIJK[1]][0];
      float offset = self->topCellOffset*abs(self->voxelDataHandler->scale[switchAxis[self->perpAxis]])*self->voxelDataHandler->voxelSize[switchAxis[self->perpAxis]]; /** add an offset to use bottom of voxel cell as height point. */

      materialExists = VoxelDataHandler_GetCurrentVoxelData( self->voxelDataHandler, (void*)&voxelData, Variable_DataType_Double );

      /** also check to see if user has defined an effective air.  If none is defined (ie, airValueMax .lt. airValueMin), this should do nothing */
      if( (voxelData > self->airValueMin || voxelData < self->airValueMax) && (self->airValueMax > self->airValueMin) )
         materialExists = False;
      if( materialExists && (*currentHeightValPtr-offset)<(float)coord[self->perpAxis] ){
         *currentHeightValPtr = (float)coord[self->perpAxis]+offset;
           *currentTestValPtr = True;
      }
   
   }
}

void _VoxelHeightFieldVariable_TestDataArray( void* voxelHeightFieldVariable ){
   VoxelHeightFieldVariable* self = (VoxelHeightFieldVariable*)voxelHeightFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   int ii, jj;
   double testHeight;

   for(ii=0; ii<self->arraySize[0]; ii++){
      for(jj=0; jj<self->arraySize[1]; jj++){
         /** take this opportunity to adjust max/min field values */
         testHeight = (double)((float***)self->voxelDataArray)[ii][jj][0];
         if( testHeight > self->localMax) self->localMax = testHeight;
         if( testHeight < self->localMin) self->localMin = testHeight;
         /** now do test */
         Journal_Firewall( ((Bool***)self->voxelDataTestArray)[ii][jj][0], self->errorStream,
            "\n\n\Error in func %s for %s '%s'\nNo height value allocate at coordinate (%f,%f).\n"
            "This suggests that no material is defined within voxel column located at this coordinate.\n"
            "Please check your voxel data, or alternatively set the MinimumHeight parameter.\n\n",
            __func__, self->type, self->name,
            self->localCrdMin[0] + ii*self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[0]]]*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[0]]],
            self->localCrdMin[1] + jj*self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[1]]]*self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[1]]]);
      }
   }
}

Bool _VoxelHeightFieldVariable_GetArrayIndexIJK( void* voxelHeightFieldVariable, Coord coord, int* indexIJK){
   VoxelHeightFieldVariable* self = (VoxelHeightFieldVariable*)voxelHeightFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   int ii;
   for(ii=0; ii<self->dim; ii++){
      indexIJK[ii] = (int) ((coord[ii] - self->localStartCrd[ii]) / ( self->voxelDataHandler->scale[switchAxis[VHFVmapAXIS[ii]]] * self->voxelDataHandler->voxelSize[switchAxis[VHFVmapAXIS[ii]]] ) ) ;
      if(indexIJK[ii]<0 || (indexIJK[ii]>=self->arraySize[ii])) return False;  /** return false to indicate that coordinate is not in local domain range */
   }
   indexIJK[2] = 0;
   return True;  /** True indicates that query location is within local domain */
}

