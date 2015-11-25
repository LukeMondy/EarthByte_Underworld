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
**  Role: Allows variable conditions to be defined on the walls of a regular mesh.
**
** Assumptions:
**
** Comments:
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Spherical_SphereBC_h__
#define __Spherical_SphereBC_h__
	

	extern const Type SphereBC_Type;
	
	extern const char* SphereBC_WallEnumToStr[SphereBC_Wall_Size];
	
	#define __SphereBC \
		/* General info */ \
		__VariableCondition \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Name						_dictionaryEntryName; \
		SphereBC_Wall				_wall; \
		Index	_entryCount; \
		WallVC_Entry*			_entryTbl; \
		Mesh*						_mesh;

	struct SphereBC { __SphereBC };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition* SphereBC_Factory(
		AbstractContext*					context,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									data );
	
	SphereBC* _SphereBC_DefaultNew( Name name );

	SphereBC* SphereBC_New(
		Name									name,
		AbstractContext*					context,
		Name									_dictionaryEntryName, 
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									_mesh );
	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define WALLVC_DEFARGS \
                VARIABLECONDITION_DEFARGS

	#define WALLVC_PASSARGS \
                VARIABLECONDITION_PASSARGS

	SphereBC* _SphereBC_New(  WALLVC_DEFARGS  );
	
	void _SphereBC_Init(
		void*	wallVC, 
		Name	_dictionaryEntryName, 
		void*	_mesh );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _SphereBC_Delete( void* wallVC );
	
	void _SphereBC_Destroy( void* wallVC, void* data );

	/* Copy */
	#define SphereBC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SphereBC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _SphereBC_Copy( void* wallVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	void _SphereBC_Build(  void* wallVC, void* data );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _SphereBC_AssignFromXML( void* wallVC, Stg_ComponentFactory* cf, void* data );
	
	void _SphereBC_BuildSelf( void* wallVC, void* data );
	
	void _SphereBC_ReadDictionary( void* variableCondition, void* dictionary );
	
	IndexSet* _SphereBC_GetSet( void* variableCondition );
	
	VariableCondition_VariableIndex _SphereBC_GetVariableCount( void* variableCondition, Index globalIndex );
	
	Variable_Index _SphereBC_GetVariableIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex );
						
	VariableCondition_ValueIndex _SphereBC_GetValueIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex	varIndex );
						
	VariableCondition_ValueIndex _SphereBC_GetValueCount( void* variableCondition );
	
	VariableCondition_Value _SphereBC_GetValue( void* variableCondition, VariableCondition_ValueIndex valIndex );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	
#endif /* __Spherical_SphereBC_h__ */

