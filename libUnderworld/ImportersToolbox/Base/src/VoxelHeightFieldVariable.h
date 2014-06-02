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

#ifndef __ImportersToolbox_Base_VoxelHeightFieldVariable_h__
#define __ImportersToolbox_Base_VoxelHeightFieldVariable_h__

	/** Textual name of this class */
	extern const Type VoxelHeightFieldVariable_Type;
	
	/** VoxelHeightFieldVariable contents */
	#define __VoxelHeightFieldVariable                                                                                         \
		__VoxelFieldVariable		                                                                                                \
		void* voxelDataTestArray;  /** array used to determine whether a value has been set at each coordinate point     */     \
		Axis perpAxis;             /** axis perpendicular to heightfield plane                                           */     \
		float minimumHeight;       /** minimum height field is initialised to. set at -HUGE_VAL by default               */     \
		double airValueMin;        /** these values (airValueMin/Max) define an 'air' material.  if the returned voxel   */     \
		double airValueMax;        /** data value is within these limits, this region is considered air, with the height */     \
		                           /** value tracing the bottom of the air region. if default, these values do nothing.  */     \
		Bool  userMinValue;        /** bool to indicate whether a user has provided a min height value                   */     \
      float topCellOffset;       /** offset of height from centre of top cell.  a value of zero results in the height being set cell centre.*/  \
                                 /**              a value of -0.5 results in the height being set at the cell bottom   */     \
                                 /**              a value of +0.5 results in the height being set at the cell top  */         \

	struct VoxelHeightFieldVariable { __VoxelHeightFieldVariable };	

	/** Creation implementation */
	void* _VoxelHeightFieldVariable_DefaultNew( Name name );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define VOXELHEIGHTFIELDVARIABLE_DEFARGS \
                VOXELFIELDVARIABLE_DEFARGS
                
	#define VOXELHEIGHTFIELDVARIABLE_PASSARGS \
                VOXELFIELDVARIABLE_PASSARGS  

	VoxelHeightFieldVariable* _VoxelHeightFieldVariable_New(  VOXELHEIGHTFIELDVARIABLE_DEFARGS  );

	/** Member initialisation implementation */
	void _VoxelHeightFieldVariable_AssignFromXML( void* VoxelHeightFieldVariable, Stg_ComponentFactory* cf, void* data ) ;
	#define _VoxelHeightFieldVariable_Build _VoxelFieldVariable_Build
	#define _VoxelHeightFieldVariable_Initialise _VoxelFieldVariable_Initialise
	#define _VoxelHeightFieldVariable_Execute _VoxelFieldVariable_Execute
	#define _VoxelHeightFieldVariable_Destroy _VoxelFieldVariable_Destroy
	
	void _VoxelHeightFieldVariable_SetArrayValue( void* voxelHeightFieldVariable, double* coord );
   void _VoxelHeightFieldVariable_SetupDataArray( void* voxelHeightFieldVariable );
   void _VoxelHeightFieldVariable_TestDataArray( void* voxelHeightFieldVariable );
   Bool _VoxelHeightFieldVariable_GetArrayIndexIJK( void* voxelHeightFieldVariable, Coord coord, int* indexIJK);
#endif /* __ImportersToolbox_Base_VoxelHeightFieldVariable_h__ */

