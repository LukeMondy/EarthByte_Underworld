/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Comments:
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Shape_LatLongRegionClass_h__
#define __StgDomain_Shape_LatLongRegionClass_h__

	/* Textual name of this class */
	extern const Type LatLongRegion_Type;

	/* LatLongRegion information */
	#define __LatLongRegion \
		/* General info */ \
		__Stg_Shape \
		/* Virtual Info */\
		\
		XYZ                                   width;            \

	struct LatLongRegion { __LatLongRegion };
	
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
   LatLongRegion* LatLongRegion_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		XYZ                                   width );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define BOX_DEFARGS \
                STG_SHAPE_DEFARGS

	#define BOX_PASSARGS \
                STG_SHAPE_PASSARGS

	LatLongRegion* _LatLongRegion_New(  BOX_DEFARGS  );
	
	void _LatLongRegion_Init( void* box, XYZ width ) ;

	/* Stg_Class_Delete LatLongRegion implementation */
	void _LatLongRegion_Delete( void* box );
	void _LatLongRegion_Print( void* box, Stream* stream );
	#define LatLongRegion_Copy( self ) \
		(LatLongRegion*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define LatLongRegion_DeepCopy( self ) \
		(LatLongRegion*) Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _LatLongRegion_Copy( void* box, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	void* _LatLongRegion_DefaultNew( Name name ) ;
	void _LatLongRegion_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) ;
	void _LatLongRegion_Build( void* box, void* data ) ;
	void _LatLongRegion_Initialise( void* box, void* data ) ;
	void _LatLongRegion_Execute( void* box, void* data );
	void _LatLongRegion_Destroy( void* box, void* data ) ;
	
	Bool _LatLongRegion_IsCoordInside( void* box, Coord coord ) ;

	double _LatLongRegion_CalculateVolume( void* box );

	void _LatLongRegion_DistanceFromCenterAxis( void* box, Coord coord, double* disVec );
	/*---------------------------------------------------------------------------------------------------------------------
	** Public member functions
	*/
	
	/*---------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/
	
	
#endif 

