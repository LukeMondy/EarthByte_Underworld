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
#include "SpatialDataFieldVariable.h"
#include "SpatialData_cppwrapper.hh"

#include <assert.h>
#include <string.h>

const Type SpatialDataFieldVariable_Type = (const Type)"SpatialDataFieldVariable";

void* _SpatialDataFieldVariable_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                      _sizeOfSelf = sizeof(SpatialDataFieldVariable);
   Type                                                              type = SpatialDataFieldVariable_Type;
   Stg_Class_DeleteFunction*                                      _delete = _FieldVariable_Delete;
   Stg_Class_PrintFunction*                                        _print = _FieldVariable_Print;
   Stg_Class_CopyFunction*                                          _copy = _FieldVariable_Copy;
   Stg_Component_DefaultConstructorFunction*          _defaultConstructor = _SpatialDataFieldVariable_DefaultNew;
   Stg_Component_ConstructFunction*                            _construct = _SpatialDataFieldVariable_AssignFromXML;
   Stg_Component_BuildFunction*                                    _build = _SpatialDataFieldVariable_Build;
   Stg_Component_InitialiseFunction*                          _initialise = _SpatialDataFieldVariable_Initialise;
   Stg_Component_ExecuteFunction*                                _execute = _FieldVariable_Execute;
   Stg_Component_DestroyFunction*                                _destroy = _SpatialDataFieldVariable_Destroy;
   AllocationType                                      nameAllocationType = NON_GLOBAL;
   FieldVariable_InterpolateValueAtFunction*          _interpolateValueAt = _SpatialDataFieldVariable_InterpolateValueAt;
   FieldVariable_GetValueFunction*            _getMinGlobalFieldMagnitude = FieldVariable_GetMinGlobalFieldMagnitude;
   FieldVariable_GetValueFunction*            _getMaxGlobalFieldMagnitude = FieldVariable_GetMaxGlobalFieldMagnitude;
   FieldVariable_CacheValuesFunction*    _cacheMinMaxGlobalFieldMagnitude = FieldVariable_CacheMinMaxGlobalFieldMagnitude;
   FieldVariable_GetCoordFunction*               _getMinAndMaxLocalCoords = FieldVariable_GetMinAndMaxLocalCoords;
   FieldVariable_GetCoordFunction*              _getMinAndMaxGlobalCoords = FieldVariable_GetMinAndMaxGlobalCoords;

#ifndef HAVE_SPATIALDATA
   Journal_Firewall( NULL, global_error_stream,
      "Error in function %s for SpatialDataFieldVariable Component:\n\n" 
      "It appears that this toolbox has not been configured with the SpatialData Library.\n"
      "Please reconfigure & recompile setting the spatialdata option:\n"
      "   ./configure.py --spatialdata-dir=/directory/where/spatialdata/build/is/\n"
      "and ensure that the configure process finds the library:\n"
      "   Checking for package spatialdata... yes", __func__, name );

#endif
   return _SpatialDataFieldVariable_New(  SPATIALDATAFIELDVARIABLE_PASSARGS  );
}

SpatialDataFieldVariable* _SpatialDataFieldVariable_New(  SPATIALDATAFIELDVARIABLE_DEFARGS  ) {
   SpatialDataFieldVariable* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(SpatialDataFieldVariable) );
   self = (SpatialDataFieldVariable*)_FieldVariable_New(  FIELDVARIABLE_PASSARGS  );

   /* Virtual functions */

   /* General info */

   /* SpatialDataFieldVariable info */
   self->filename   = NULL;
   self->approxType = NULL;
   self->fieldName  = NULL;

   return self;
}

void _SpatialDataFieldVariable_AssignFromXML( void* _SpatialDataFieldVariable, Stg_ComponentFactory* cf, void* data ) {
   SpatialDataFieldVariable* self = (SpatialDataFieldVariable*)_SpatialDataFieldVariable;

   _FieldVariable_AssignFromXML( self, cf, data );

   /** The name of the file containing the voxel data */
   self->filename   = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"Filename", NULL ));
   /** type of approximation to use.  accepts 'nearest' or 'linear' (actually, defaults to nearest if not 'linear') */
   self->approxType = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"ApproximationType", "nearest" ));
   /** field name as found in the spatialDB file.  only a single field currently allowed */
   self->fieldName  = StG_Strdup(Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"DBFieldName", NULL ));

   Journal_Firewall( self->filename, global_error_stream,
      "Error in function %s for SpatialDataFieldVariable Component: No filename provided." 
      "\nPlease set the parameter 'Filename' for the '%s' component.", __func__, self->name );
   Journal_Firewall( self->fieldName, global_error_stream,
      "Error in function %s for SpatialDataFieldVariable Component: No SpatialDB field name provided." 
      "\nPlease set the parameter 'DBFieldName' for the '%s' component.", __func__, self->name );

}

void _SpatialDataFieldVariable_Build( void* _SpatialDataFieldVariable, void* data ) {
   SpatialDataFieldVariable* self = (SpatialDataFieldVariable*)_SpatialDataFieldVariable;
   
   _FieldVariable_Build( self, data );
#ifdef HAVE_SPATIALDATA
   _SpatialData_cppwrapper_Build( &(self->cppdata), self->filename, self->approxType, self->fieldName );
#endif

}

void _SpatialDataFieldVariable_Initialise( void* _SpatialDataFieldVariable, void* data ) {
   SpatialDataFieldVariable* self = (SpatialDataFieldVariable*)_SpatialDataFieldVariable;
   
   _FieldVariable_Initialise( self, data );
   
}

void _SpatialDataFieldVariable_Execute( void* SpatialDataFieldVariable, void* data ) {

}

void _SpatialDataFieldVariable_Destroy( void* _SpatialDataFieldVariable, void* data ) {
   SpatialDataFieldVariable* self = (SpatialDataFieldVariable*)_SpatialDataFieldVariable;

#ifdef HAVE_SPATIALDATA
   _SpatialData_cppwrapper_Destroy( &(self->cppdata) );
#endif
   Memory_Free(self->filename);
   Memory_Free(self->approxType);
   Memory_Free(self->fieldName);
}

InterpolationResult _SpatialDataFieldVariable_InterpolateValueAt( void* _SpatialDataFieldVariable, double* coord, double* value ){
   SpatialDataFieldVariable* self = (SpatialDataFieldVariable*)_SpatialDataFieldVariable;   
   int returned;
#ifdef HAVE_SPATIALDATA
   returned = _SpatialData_cppwrapper_evaluate( &(self->cppdata), coord, value );
#endif
   if (returned == 0)
      return LOCAL;
   else 
      return OUTSIDE_GLOBAL;
}

