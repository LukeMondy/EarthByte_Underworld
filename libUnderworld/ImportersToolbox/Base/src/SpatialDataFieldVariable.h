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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __ImportersToolbox_Base_SpatialDataFieldVariable_h__
#define __ImportersToolbox_Base_SpatialDataFieldVariable_h__

	/** Textual name of this class */
	extern const Type SpatialDataFieldVariable_Type;
	
	/** SpatialDataFieldVariable contents */
	#define __SpatialDataFieldVariable                          \
		/* General info */                                      \
		__FieldVariable                                         \
		                                                        \
		/* Member info */                                       \
		Name                     filename;                      \
		Name                     approxType;                    \
		Name                     fieldName;                     \
        void*                    cppdata;

	struct SpatialDataFieldVariable { __SpatialDataFieldVariable };	

	/** Creation implementation */
	void* _SpatialDataFieldVariable_DefaultNew( Name name );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

   #define SPATIALDATAFIELDVARIABLE_DEFARGS \
                FIELDVARIABLE_DEFARGS
                
   #define SPATIALDATAFIELDVARIABLE_PASSARGS \
                FIELDVARIABLE_PASSARGS  

   SpatialDataFieldVariable* _SpatialDataFieldVariable_New(  SPATIALDATAFIELDVARIABLE_DEFARGS  );

	void _SpatialDataFieldVariable_AssignFromXML( void* SpatialDataFieldVariable, Stg_ComponentFactory* cf, void* data ) ;

	void _SpatialDataFieldVariable_Build( void* SpatialDataFieldVariable, void* data ) ;

	void _SpatialDataFieldVariable_Execute( void* SpatialDataFieldVariable, void* data ) ;

	void _SpatialDataFieldVariable_Destroy( void* SpatialDataFieldVariable, void* data ) ;

	void _SpatialDataFieldVariable_Initialise( void* SpatialDataFieldVariable, void* data ) ;
	
	InterpolationResult _SpatialDataFieldVariable_InterpolateValueAt( void* SpatialDataFieldVariable, double* coord, double* value );

#endif /* __ImportersToolbox_Base_SpatialDataFieldVariable_h__ */

