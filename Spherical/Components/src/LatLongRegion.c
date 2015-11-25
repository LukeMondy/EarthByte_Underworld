/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) )
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "LatLongRegion.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type LatLongRegion_Type = "LatLongRegion";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

LatLongRegion* _LatLongRegion_New(  BOX_DEFARGS  )
{
	LatLongRegion* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(LatLongRegion) );
	self = (LatLongRegion*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */

	return self;
}

void _LatLongRegion_Init( void* shape, XYZ width ) {
	LatLongRegion* self = (LatLongRegion*)shape;
	
	memcpy( self->width, width, sizeof(XYZ));
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _LatLongRegion_Delete( void* shape ) {
	LatLongRegion* self = (LatLongRegion*)shape;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}


void _LatLongRegion_Print( void* shape, Stream* stream ) {
	LatLongRegion* self = (LatLongRegion*)shape;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _LatLongRegion_Copy( void* shape, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	LatLongRegion*	self = (LatLongRegion*)shape;
	LatLongRegion*	newLatLongRegion;
	
	newLatLongRegion = (LatLongRegion*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );
	
	memcpy( newLatLongRegion->width, self->width, sizeof(XYZ));
	
	return (void*)newLatLongRegion;
}

void* _LatLongRegion_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(LatLongRegion);
	Type                                                          type = LatLongRegion_Type;
	Stg_Class_DeleteFunction*                                  _delete = _LatLongRegion_Delete;
	Stg_Class_PrintFunction*                                    _print = _LatLongRegion_Print;
	Stg_Class_CopyFunction*                                      _copy = _LatLongRegion_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _LatLongRegion_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _LatLongRegion_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _LatLongRegion_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _LatLongRegion_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _LatLongRegion_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _LatLongRegion_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _LatLongRegion_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _LatLongRegion_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _LatLongRegion_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _LatLongRegion_New(  BOX_PASSARGS  );
}


void _LatLongRegion_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) {
	LatLongRegion*	             self          = (LatLongRegion*) shape;
	Dictionary*          dictionary    = Dictionary_GetDictionary( cf->componentDict, self->name );
	XYZ                  width;
	double               start, end;
	Dictionary_Entry_Key startKey      = StG_Strdup("start");
	Dictionary_Entry_Key endKey        = StG_Strdup("end");
	Dimension_Index      dim_I;

        char* words[]={"Radius", "Long", "Lat"};
        char destStart[500];
        char destEnd[500];


	_Stg_Shape_AssignFromXML( self, cf, data );

        for( dim_I=0; dim_I<3; dim_I++ ) {
           // make start and end dictionary words
           strcpy( destStart, startKey ); strcat( destStart, words[dim_I] );
           strcpy( destEnd, endKey ); strcat( destEnd, words[dim_I] );

           // get from the dictionary
           Stg_ComponentFactory_GetRequiredDouble( cf, self->name, (Dictionary_Entry_Key) destStart, &start ); 
           Stg_ComponentFactory_GetRequiredDouble( cf, self->name, (Dictionary_Entry_Key) destEnd, &end ); 

           width[dim_I] = fabs(end - start);
           if( end > start ) self->centre[ dim_I ] = start + 0.5 * width[dim_I];
           else self->centre[ dim_I ] = start - 0.5 * width[dim_I];
        }

	Memory_Free( startKey );
	Memory_Free( endKey );

	_LatLongRegion_Init( self, width );
}

void _LatLongRegion_Build( void* shape, void* data ) {
	LatLongRegion*	self = (LatLongRegion*)shape;

	_Stg_Shape_Build( self, data );
}
void _LatLongRegion_Initialise( void* shape, void* data ) {
	LatLongRegion*	self = (LatLongRegion*)shape;
	
	_Stg_Shape_Initialise( self, data );
}
void _LatLongRegion_Execute( void* shape, void* data ) {
	LatLongRegion*	self = (LatLongRegion*)shape;
	
	_Stg_Shape_Execute( self, data );
}
void _LatLongRegion_Destroy( void* shape, void* data ) {
	LatLongRegion*	self = (LatLongRegion*)shape;
    
	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _LatLongRegion_IsCoordInside( void* shape, Coord coord ) {
   LatLongRegion* self = (LatLongRegion*)shape;
   Coord newCoord;
   Dimension_Index dim_I;
   Coord rtp;

   if( self->dim == 2 ) {
     Spherical_XYZ2RTP2D( coord, rtp );
     rtp[1] = rtp[1] * (180.0/M_PI); // convert to degrees
   }
   else {
     /* Transform coordinate from x,y,z into radius,theta,phi */
     Spherical_XYZ2RTP3D( coord, rtp );
     rtp[1] = rtp[1] * (180.0/M_PI); // convert to degrees
     rtp[2] = rtp[2] * (180.0/M_PI); // convert to degrees
   }

   for( dim_I=0; dim_I<self->dim; dim_I++ ) {
      // offset rtp coords into shape reference frame 
      rtp[dim_I] = rtp[dim_I] - self->centre[dim_I];
      // test if new rtp point is within rtp widths
      if( fabs(rtp[dim_I]) > 0.5*self->width[dim_I] ) {
         return False;
      }
   }

   return True;
}

double _LatLongRegion_CalculateVolume( void* shape ) {
   Journal_Firewall( False, global_error_stream,
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}

void _LatLongRegion_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ) {
	/* To be implemented */
	Journal_Firewall( False, global_error_stream,
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}



