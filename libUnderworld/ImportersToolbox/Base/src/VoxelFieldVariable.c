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
#include "VoxelDataHandler_ndarray.h"
#include "VoxelFieldVariable.h"

#include <assert.h>
#include <string.h>

const Type VoxelFieldVariable_Type = "VoxelFieldVariable";

void* _VoxelFieldVariable_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                         _sizeOfSelf = sizeof(VoxelFieldVariable);
   Type                                                                 type = VoxelFieldVariable_Type;
   Stg_Class_DeleteFunction*                                         _delete = _FieldVariable_Delete;
   Stg_Class_PrintFunction*                                           _print = _FieldVariable_Print;
   Stg_Class_CopyFunction*                                             _copy = _FieldVariable_Copy;
   Stg_Component_DefaultConstructorFunction*             _defaultConstructor = _VoxelFieldVariable_DefaultNew;
   Stg_Component_ConstructFunction*                               _construct = _VoxelFieldVariable_AssignFromXML;
   Stg_Component_BuildFunction*                                       _build = _VoxelFieldVariable_Build;
   Stg_Component_InitialiseFunction*                             _initialise = _VoxelFieldVariable_Initialise;
   Stg_Component_ExecuteFunction*                                   _execute = _VoxelFieldVariable_Execute;
   Stg_Component_DestroyFunction*                                   _destroy = _VoxelFieldVariable_Destroy;
   AllocationType                                         nameAllocationType = NON_GLOBAL;
   FieldVariable_InterpolateValueAtFunction*             _interpolateValueAt = _VoxelFieldVariable_InterpolateValueAt;
   FieldVariable_GetValueFunction*               _getMinGlobalFieldMagnitude = _VoxelFieldVariable_GetMinGlobalFieldMagnitude;
   FieldVariable_GetValueFunction*               _getMaxGlobalFieldMagnitude = _VoxelFieldVariable_GetMaxGlobalFieldMagnitude;
   FieldVariable_CacheValuesFunction*       _cacheMinMaxGlobalFieldMagnitude = _VoxelFieldVariable_CacheMinMaxGlobalFieldMagnitude;
   FieldVariable_GetCoordFunction*                  _getMinAndMaxLocalCoords = _VoxelFieldVariable_GetMinAndMaxLocalCoords;
   FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords = _VoxelFieldVariable_GetMinAndMaxGlobalCoords;
   VoxelFieldVariable_SetupDataArray*                        _setupDataArray = _VoxelFieldVariable_SetupDataArray;
   VoxelFieldVariable_TestDataArray*                          _testDataArray = _VoxelFieldVariable_TestDataArray;
   VoxelFieldVariable_SetArrayValue*                          _setArrayValue = _VoxelFieldVariable_SetArrayValue;
   VoxelFieldVariable_GetArrayIndexIJK*                    _getArrayIndexIJK = _VoxelFieldVariable_GetArrayIndexIJK;
   
   return _VoxelFieldVariable_New(  VOXELFIELDVARIABLE_PASSARGS  );
}

VoxelFieldVariable* _VoxelFieldVariable_New(  VOXELFIELDVARIABLE_DEFARGS  ) {
   VoxelFieldVariable* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(VoxelFieldVariable) );
   self = (VoxelFieldVariable*)_FieldVariable_New(  FIELDVARIABLE_PASSARGS  );

   self->_setupDataArray   = _setupDataArray;
   self->_testDataArray    = _testDataArray;
   self->_setArrayValue    = _setArrayValue;
   self->_getArrayIndexIJK = _getArrayIndexIJK;
   
   self->localMax = -HUGE_VAL;
   self->localMin =  HUGE_VAL;

   self->decompositionMesh = NULL;
   self->voxelDataArray    = NULL;

   self->minMaxCached    = False;

   return self;
}

void _VoxelFieldVariable_AssignFromXML( void* _VoxelFieldVariable, Stg_ComponentFactory* cf, void* data ) {
   VoxelFieldVariable*         self = (VoxelFieldVariable*)_VoxelFieldVariable;

   _FieldVariable_AssignFromXML( self, cf, data );

   /** voxel data handler provides all the required information about the voxel dataset */
   self->voxelDataHandler = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"VoxelDataHandler", VoxelDataHandler_Abstract, True, data  );
   /** the feMesh is used for parallel computation, so that only the required parts of the voxel dataset are loaded for each processor. */
   /** if no decomposition mesh is provided, the entire dataset is loaded on each processor */
   self->decompositionMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"decompositionMesh", Mesh, False, data  );
   /** for 2d simulations, select which Z slice to take from voxel data (which are natively 3d).  defaults to zero */
   self->zSliceCoord = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"zSliceCoord", 0 );
   /** if querying outside the domain, should we used nearest cell? */
   self->useNearestCell = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"UseNearestCellIfOutside", False );

   self->errorStream       = Journal_Register( Error_Type, (Name)self->type  );
   self->dataType          = self->voxelDataHandler->dataType;

   /* if using ndarray, need to set this flag as memory is allocated within numpy */
   if (Stg_Class_CompareType( self->voxelDataHandler, VoxelDataHandler_ndarray_Type)){
      self->useExternalArray = True;
   }

   self->fieldComponentCount = 1;

}

void _VoxelFieldVariable_Build( void* _VoxelFieldVariable, void* data ) {
   VoxelFieldVariable* self = (VoxelFieldVariable*)_VoxelFieldVariable;

   /** build parent */
   _FieldVariable_Build( _VoxelFieldVariable, data );

   Stg_Component_Build( self->decompositionMesh, data, False );
   Stg_Component_Build( self->voxelDataHandler, data, False );
}


void _VoxelFieldVariable_Initialise( void* _VoxelFieldVariable, void* data ) {
   VoxelFieldVariable* self = (VoxelFieldVariable*)_VoxelFieldVariable;
   int voxelCount;

   /** initialise parent */
   _FieldVariable_Initialise( _VoxelFieldVariable, data );
   
   Stg_Component_Initialise( self->voxelDataHandler, data, False );
   Stg_Component_Initialise( self->decompositionMesh, data, False );
   
   /** now that all required components are initialised, we proceed to build/init our voxel field variable. */
   /** first determine and setup array for attributes storage */
   self->_setupDataArray( _VoxelFieldVariable );

   /** now initialise the data */
   /** get voxel count */
   VoxelDataHandler_GetTotalVoxelCount( self->voxelDataHandler, &voxelCount );
   /** initialise voxel data handler */
   VoxelDataHandler_GotoVoxelDataTop( self->voxelDataHandler );

   if( !self->useExternalArray ){
      int vox_i;
      Progress* prog;
      char* titleString;

      /* Use a progress meter. */
      prog = Progress_New();
      Progress_SetStream( prog, Journal_MyStream( Info_Type, self ) );
      Stg_asprintf( &titleString, "%s - Parsing Voxel Data File (%s)", self->type, self->voxelDataHandler->filename );
      Progress_SetTitle( prog, titleString );
      Progress_SetRange( prog, 0, voxelCount );

      for(vox_i=0; vox_i<voxelCount; vox_i++){
         Coord coord;
         VoxelDataHandler_GetCurrentVoxelCoord( (void*)self->voxelDataHandler, (double*)&coord );
         
         self->_setArrayValue( (void*)self, (double*)&coord );

         /** increment voxel data cursor */
         Progress_Increment( prog );
         VoxelDataHandler_IncrementVoxelIndex( self->voxelDataHandler );
      }
      /** if required, perform any tests to ensure correct initialisation of data */
      self->_testDataArray( _VoxelFieldVariable );

      /** force re-cache of min and max */
      self->minMaxCached = False;
      _VoxelFieldVariable_GetMinGlobalFieldMagnitude( self );

      Memory_Free( titleString );
      /* Delete progress meter. */
      Stg_Class_Delete( prog );
   }

}

void _VoxelFieldVariable_Destroy( void* _VoxelFieldVariable, void* data ) {
   VoxelFieldVariable* self = (VoxelFieldVariable*)_VoxelFieldVariable;

   /** destroy parent */
   _FieldVariable_Destroy( _VoxelFieldVariable, data );

   Stg_Component_Destroy( self->decompositionMesh, data, False );
   Stg_Component_Destroy( self->voxelDataHandler, data, False );

}

double _VoxelFieldVariable_GetMinGlobalFieldMagnitude( void* voxelFieldVariable ) {
   return FieldVariable_GetMinGlobalFieldMagnitude(voxelFieldVariable);
}

double _VoxelFieldVariable_GetMaxGlobalFieldMagnitude( void* voxelFieldVariable ) {
   return FieldVariable_GetMaxGlobalFieldMagnitude(voxelFieldVariable);
}

void _VoxelFieldVariable_CacheMinMaxGlobalFieldMagnitude( void* voxelFieldVariable ) {
   VoxelFieldVariable* self = (VoxelFieldVariable*) voxelFieldVariable;
   if(self->minMaxCached == False){
      int ii;
      for(ii=0; ii<self->arraySize[0]*self->arraySize[1]*self->arraySize[2]; ii++) {
         double value;
         switch ( self->dataType )
         {
            case Variable_DataType_Char   :  value = (double) *(signed char*)(self->voxelDataArray + ii*sizeof(signed char));   break;
            case Variable_DataType_Int    :  value = (double) *(        int*)(self->voxelDataArray + ii*sizeof(        int));   break;
            case Variable_DataType_Float  :  value = (double) *(      float*)(self->voxelDataArray + ii*sizeof(      float));   break;
            case Variable_DataType_Double :  value = (double) *(     double*)(self->voxelDataArray + ii*sizeof(     double));   break;
         }

         if( value > self->localMax) self->localMax = value;
         if( value < self->localMin) self->localMin = value;

      }

      MPI_Comm comm = MPI_COMM_WORLD;
      if (self->context) comm = self->context->communicator;
      /** Find upper and lower bounds on all processors */
      (void)MPI_Allreduce( &self->localMin, &self->magnitudeMin, 1, MPI_DOUBLE, MPI_MIN, comm);
      (void)MPI_Allreduce( &self->localMax, &self->magnitudeMax, 1, MPI_DOUBLE, MPI_MAX, comm);
      self->minMaxCached = True;
   }
}

void _VoxelFieldVariable_GetMinAndMaxLocalCoords( void* voxelFieldVariable, double* min, double* max ) {
   VoxelFieldVariable*  self = (VoxelFieldVariable*)voxelFieldVariable;
   int ii;

   for(ii=0; ii<self->dim; ii++){
      min[ii] = self->localCrdMin[ii];
      max[ii] = self->localCrdMax[ii];
   }
}

void _VoxelFieldVariable_GetMinAndMaxGlobalCoords( void* voxelFieldVariable, double* min, double* max ) {
   VoxelFieldVariable*  self = (VoxelFieldVariable*)voxelFieldVariable;
   MPI_Comm comm = MPI_COMM_WORLD;
   if (self->context) comm = self->context->communicator;

   /** Find upper and lower bounds on all processors */
   int ii;
   for(ii=0; ii<self->dim; ii++){
      (void)MPI_Allreduce( &self->localCrdMin[ii], &min[ii], 1, MPI_DOUBLE, MPI_MIN, comm);
      (void)MPI_Allreduce( &self->localCrdMax[ii], &max[ii], 1, MPI_DOUBLE, MPI_MAX, comm);
   }
}

InterpolationResult _VoxelFieldVariable_InterpolateValueAt( void* voxelFieldVariable, double* globalCoord, double* value ) {
   VoxelFieldVariable* self = (VoxelFieldVariable*)voxelFieldVariable;
   InterpolationResult  retValue = OUTSIDE_GLOBAL;
   Bool                 local;
   int                  indexIJK[3] = { 0, 0, 0 };
   Coord                coord;

   coord[0] = globalCoord[0];
   coord[1] = globalCoord[1];
   if(self->dim == 2)
      coord[2] = self->zSliceCoord;
   else
      coord[2] = globalCoord[2];

   local = self->_getArrayIndexIJK( (void*)self, coord, (int*)&indexIJK);

   if(local){
      int posy = indexIJK[0] + self->arraySize[0]*indexIJK[1] + self->arraySize[0]*self->arraySize[1]*indexIJK[2];

      switch ( self->dataType )
      {
         case Variable_DataType_Char   :  *value = (double) *(signed char*)( self->voxelDataArray + posy*sizeof(signed char) );   break;
         case Variable_DataType_Int    :  *value = (double) *(        int*)( self->voxelDataArray + posy*sizeof(        int) );   break;
         case Variable_DataType_Float  :  *value = (double) *(      float*)( self->voxelDataArray + posy*sizeof(      float) );   break;
         case Variable_DataType_Double :  *value = (double) *(     double*)( self->voxelDataArray + posy*sizeof(     double) );   break;
      }
      retValue = LOCAL;
   }

   return retValue;
}

void _VoxelFieldVariable_SetupDataArray( void* voxelFieldVariable ){
   VoxelFieldVariable* self = (VoxelFieldVariable*)voxelFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   if( self->decompositionMesh ){
      Journal_Firewall( NULL, self->errorStream, "\n\nError in %s for %s '%s' - While everything is in place, using a decomposition mesh has not been tested! Use with caution.\n\n", __func__, self->type, self->name );
      double  displacement[3];
      int     ii;
      double meshCrdMin[3], meshCrdMax[3];

      /** first ensure all information is updated if a deformation has occcured */
      Mesh_DeformationUpdate( self->decompositionMesh );
      Mesh_GetLocalCoordRange( self->decompositionMesh, &meshCrdMin, &meshCrdMax );
      for(ii=0; ii<3; ii++){
         /** determine range of voxels required over this processors domain */
         self->arraySize[ii] = (int) ( meshCrdMax[ii] - meshCrdMin[ii] ) / ( abs(self->voxelDataHandler->scale[switchAxis[ii]])*self->voxelDataHandler->voxelSize[switchAxis[ii]]) + 1;
         /** get displacement of mesh from start of voxel data depending on which way we are counting */
         if(self->voxelDataHandler->scale[switchAxis[ii]] > 0){
            displacement[ii] = meshCrdMin[ii] - self->voxelDataHandler->startCrd[switchAxis[ii]] - 0.5*self->voxelDataHandler->scale[switchAxis[ii]]*self->voxelDataHandler->voxelSize[switchAxis[ii]];
            /** determine local support coords */
            self->localCrdMin[ii]   = self->voxelDataHandler->minDomainCoordsXYZ[ii] + (int) displacement[ii];
            self->localCrdMax[ii]   = self->localCrdMin[ii] + self->arraySize[ii]*abs(self->voxelDataHandler->scale[switchAxis[ii]])*self->voxelDataHandler->voxelSize[switchAxis[ii]];
            self->localStartCrd[ii] = self->localCrdMin[ii]; 
         } else {
            displacement[ii] = self->voxelDataHandler->startCrd[switchAxis[ii]] + 0.5*self->voxelDataHandler->scale[switchAxis[ii]]*self->voxelDataHandler->voxelSize[switchAxis[ii]] - meshCrdMax[ii];
            /** determine local support coords */
            self->localCrdMax[ii]    = self->voxelDataHandler->maxDomainCoordsXYZ[ii] - (int) displacement[ii];
            self->localCrdMin[ii]    = self->localCrdMax[ii] - self->arraySize[ii]*abs(self->voxelDataHandler->scale[switchAxis[ii]])*self->voxelDataHandler->voxelSize[switchAxis[ii]];
            self->localStartCrd[ii]  = self->localCrdMax[ii];
         }
         Journal_Firewall( displacement[ii] > 0, self->errorStream, "\n\nError in %s for %s '%s' - Voxel data is defined over a smaller region than local domain, which is not valid.\n\n", __func__, self->type, self->name );
      }
   } else {
      int ii;
      for(ii=0; ii<3; ii++){
         /** local array size equal to the voxel size */
         self->arraySize[ii] = self->voxelDataHandler->numCells[switchAxis[ii]];
         /** local support will be set to the support of the voxel dataset */
         self->localCrdMin[ii] = self->voxelDataHandler->minDomainCoordsXYZ[ii];
         self->localCrdMax[ii] = self->voxelDataHandler->maxDomainCoordsXYZ[ii];
         /** check for counting direction */
         if(self->voxelDataHandler->scale[switchAxis[ii]] > 0){
            self->localStartCrd[ii] = self->localCrdMin[ii]; 
         } else {
            self->localStartCrd[ii] = self->localCrdMax[ii];
         }
      }
   }

   if (!self->useExternalArray) {
      if( self->voxelDataArray != NULL )
         Memory_Free( self->voxelDataArray );

      int totCells=self->arraySize[0]*self->arraySize[1]*self->arraySize[2];
      /** allocate memory */
      switch ( self->dataType )
      {
         case Variable_DataType_Char   : self->voxelDataArray = Memory_Alloc_Bytes_Unnamed( sizeof(signed char)*totCells, (Name)"voxelfield array" );   break;
         case Variable_DataType_Int    : self->voxelDataArray = Memory_Alloc_Bytes_Unnamed( sizeof(        int)*totCells, (Name)"voxelfield array" );   break;
         case Variable_DataType_Float  : self->voxelDataArray = Memory_Alloc_Bytes_Unnamed( sizeof(      float)*totCells, (Name)"voxelfield array" );   break;
         case Variable_DataType_Double : self->voxelDataArray = Memory_Alloc_Bytes_Unnamed( sizeof(     double)*totCells, (Name)"voxelfield array" );   break;
         default                       : Journal_Firewall( NULL, self->errorStream, "\n\nError in %s for %s '%s' - Datatype not know or not supported.\n\n",__func__,self->type,self->name );  break;
      }
   } else {
      self->voxelDataArray = (void*)((VoxelDataHandler_ndarray*)self->voxelDataHandler)->ndPointer;
   }
   
}

void _VoxelFieldVariable_TestDataArray( void* voxelFieldVariable ){};

void _VoxelFieldVariable_SetArrayValue( void* voxelFieldVariable, double* coord ){
   VoxelFieldVariable* self = (VoxelFieldVariable*)voxelFieldVariable;
   int indexIJK[3] = { 0, 0, 0 };
   
   /** if coord is in local range, set value */
   if( self->_getArrayIndexIJK( voxelFieldVariable, coord, &indexIJK) ){
      double doubleData;
      /** allocate memory.. simply choose size double to ensure enough memory */
      void* voxelData = Memory_Alloc_Bytes_Unnamed( sizeof(double), "double" );
      VoxelDataHandler_GetCurrentVoxelData( self->voxelDataHandler, voxelData, self->dataType );

      int posy = indexIJK[0] + self->arraySize[0]*indexIJK[1] + self->arraySize[0]*self->arraySize[1]*indexIJK[2];

      switch ( self->dataType )
      {
         case Variable_DataType_Char   :  *(( signed char*)( self->voxelDataArray + posy*sizeof(signed char) )) = *((signed char*)voxelData);  doubleData = (double)*((signed char*)voxelData);  break;
         case Variable_DataType_Int    :  *((         int*)( self->voxelDataArray + posy*sizeof(        int) )) =         *((int*)voxelData);  doubleData = (double)*((        int*)voxelData);  break;
         case Variable_DataType_Float  :  *((       float*)( self->voxelDataArray + posy*sizeof(      float) )) =       *((float*)voxelData);  doubleData = (double)*((      float*)voxelData);  break;
         case Variable_DataType_Double :  *((      double*)( self->voxelDataArray + posy*sizeof(     double) )) =      *((double*)voxelData);  doubleData = (double)*((     double*)voxelData);  break;
      }

      if( doubleData > self->localMax) self->localMax = doubleData;
      if( doubleData < self->localMin) self->localMin = doubleData;
      
      Memory_Free( voxelData );
   }
}

Bool _VoxelFieldVariable_GetArrayIndexIJK( void* voxelFieldVariable, Coord coord, int* indexIJK){
   VoxelFieldVariable* self = (VoxelFieldVariable*)voxelFieldVariable;
   Index* switchAxis = &self->voxelDataHandler->switchAxis[0];
   int ii;
   for(ii=0; ii<self->dim; ii++){
      indexIJK[ii] = (int) ((coord[ii] - self->localStartCrd[ii])/ ( self->voxelDataHandler->scale[switchAxis[ii]]*self->voxelDataHandler->voxelSize[switchAxis[ii]] ) );
      if(indexIJK[ii]<0 || (indexIJK[ii]>=self->arraySize[ii])){
         if ( !self->useNearestCell ) {
            return False;  /** return false to indicate that coordinate is not in local domain range */
         } else {
            if( indexIJK[ii]<0 )
               indexIJK[ii] = 0;   /** ok, we are less than min, so use min instead */
            else
               indexIJK[ii] = self->arraySize[ii]-1;  /** ok, we must be more than max, so use max instead */
         }
      }
   }
   return True;  /** True indicates that query location is within local domain */
}

