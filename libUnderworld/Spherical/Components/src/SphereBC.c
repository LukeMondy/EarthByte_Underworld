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
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/StgDomain.h>

#include "Components.h"

#include <string.h>
#include <assert.h>


const Type SphereBC_Type = "SphereBC";
const Name defaultSphereBCName = "defaultSphereBCName";

const char* SphereBC_WallEnumToStr[SphereBC_Wall_Size] = {
    "Inner",
    "Outer"
};


/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* SphereBC_Factory(
    AbstractContext*					context,
    Variable_Register*				variable_Register,
    ConditionFunction_Register*	conFunc_Register,
    Dictionary*							dictionary,
    void*									data )
{
    return (VariableCondition*)SphereBC_New( defaultSphereBCName, context, NULL, variable_Register, conFunc_Register, dictionary, (Mesh*)data );
}

SphereBC* SphereBC_New(
    Name									name,
    AbstractContext*					context,
    Name									_dictionaryEntryName,
    Variable_Register*				variable_Register,
    ConditionFunction_Register*	conFunc_Register,
    Dictionary*							dictionary,
    void*									_mesh )
{
    SphereBC* self = _SphereBC_DefaultNew( name );

    self->isConstructed = True;
    _VariableCondition_Init( self, context, variable_Register, conFunc_Register, dictionary );
    _SphereBC_Init( self, _dictionaryEntryName, _mesh );

    return self;
}

SphereBC* _SphereBC_DefaultNew( Name name )
{
    /* Variables set in this function */
    SizeT                                               _sizeOfSelf = sizeof(SphereBC);
    Type                                                       type = SphereBC_Type;
    Stg_Class_DeleteFunction*                               _delete = _SphereBC_Delete;
    Stg_Class_PrintFunction*                                 _print = NULL;
    Stg_Class_CopyFunction*                                   _copy = _SphereBC_Copy;
    Stg_Component_DefaultConstructorFunction*   _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)_SphereBC_DefaultNew;
    Stg_Component_ConstructFunction*                     _construct = _SphereBC_AssignFromXML;
    Stg_Component_BuildFunction*                             _build = _SphereBC_Build;
    Stg_Component_InitialiseFunction*                   _initialise = _VariableCondition_Initialise;
    Stg_Component_ExecuteFunction*                         _execute = _VariableCondition_Execute;
    Stg_Component_DestroyFunction*                         _destroy = _SphereBC_Destroy;
    AllocationType                               nameAllocationType = NON_GLOBAL;
    VariableCondition_BuildSelfFunc*                     _buildSelf = _SphereBC_BuildSelf;
    VariableCondition_PrintConciseFunc*               _printConcise = NULL;
    VariableCondition_ReadDictionaryFunc*           _readDictionary = _SphereBC_ReadDictionary;
    VariableCondition_GetSetFunc*                           _getSet = _SphereBC_GetSet;
    VariableCondition_GetVariableCountFunc*       _getVariableCount = _SphereBC_GetVariableCount;
    VariableCondition_GetVariableIndexFunc*       _getVariableIndex = _SphereBC_GetVariableIndex;
    VariableCondition_GetValueIndexFunc*             _getValueIndex = _SphereBC_GetValueIndex;
    VariableCondition_GetValueCountFunc*             _getValueCount = _SphereBC_GetValueCount;
    VariableCondition_GetValueFunc*                       _getValue = _SphereBC_GetValue;
    VariableCondition_ApplyFunc*                             _apply = _VariableCondition_Apply;

    return _SphereBC_New(  WALLVC_PASSARGS  );
}

SphereBC* _SphereBC_New(  WALLVC_DEFARGS  )
{
    SphereBC* self;

    /* Allocate memory/General info */
    assert( _sizeOfSelf >= sizeof(SphereBC) );
    self = (SphereBC*)_VariableCondition_New(  VARIABLECONDITION_PASSARGS  );

    /* Virtual info */

    /* Stg_Class info */

    return self;
}

void _SphereBC_Init( void* wallVC, Name _dictionaryEntryName, void* _mesh )
{
    SphereBC* self = (SphereBC*)wallVC;

    self->_dictionaryEntryName = _dictionaryEntryName;
    self->_mesh = (Mesh*)_mesh;
    self->_wall = SphereBC_Wall_Size;
    self->_entryTbl = 0;
    self->_entryCount = 0;
}


/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _SphereBC_ReadDictionary( void* variableCondition, void* dictionary )
{
    SphereBC*			self = (SphereBC*)variableCondition;
    Dictionary_Entry_Value*	vcDictVal;
    Dictionary_Entry_Value	_vcDictVal;
    Dictionary_Entry_Value*	varsVal;
    WallVC_Entry_Index	entry_I;

    /* Find dictionary entry */
    if (self->_dictionaryEntryName)
        vcDictVal = Dictionary_Get( dictionary, (Dictionary_Entry_Key)self->_dictionaryEntryName );
    else {
        vcDictVal = &_vcDictVal;
        Dictionary_Entry_Value_InitFromStruct(vcDictVal, dictionary);
    }

    if (vcDictVal) {
        char*	wallStr;

        /* Obtain which wall */
        wallStr = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"wall" ) );
        if (     !strcasecmp(wallStr, "inner") )
            self->_wall = SphereBC_Wall_Inner;
        else if (!strcasecmp(wallStr, "outer") )
            self->_wall = SphereBC_Wall_Outer;
        else {
            Stream*	errorStr = Journal_Register( Error_Type, (Name)self->type  );
            Journal_Firewall( 0 , errorStr, "Error- in %s: wallVC type '%s' is invalid.  Valid wallVC types are: \n\n"
                              "   MinK (back)\n   MinI (left)\n   MinJ (bottom)\n   MaxI (right)\n   MaxJ (top)\n   MaxK (front)\n   MinIMinJ (bottomLeft)\n   MaxIMinJ (bottomRight)\n\n"
                              , __func__, wallStr );

            self->_wall = SphereBC_Wall_Size; /* invalid entry */
        }

        /* Obtain the variable entries */
        self->_entryCount = Dictionary_Entry_Value_GetCount(Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"variables") );
        self->_entryTbl = Memory_Alloc_Array( WallVC_Entry, self->_entryCount, "SphereBC->_entryTbl" );
        varsVal = Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"variables");

        for (entry_I = 0; entry_I < self->_entryCount; entry_I++ ) {
            char*			valType;
            Dictionary_Entry_Value*	valueEntry;
            Dictionary_Entry_Value*	varDictListVal;

            varDictListVal = Dictionary_Entry_Value_GetElement(varsVal, entry_I);
            valueEntry = Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"value" );

            self->_entryTbl[entry_I].varName = Dictionary_Entry_Value_AsString(
                                                   Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"name") );

            valType = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"type") );
            if (0 == strcasecmp(valType, "func")) {
                char*	funcName = Dictionary_Entry_Value_AsString(valueEntry);
                Index	cfIndex;

                self->_entryTbl[entry_I].value.type = VC_ValueType_CFIndex;
                cfIndex = ConditionFunction_Register_GetIndex( self->conFunc_Register, funcName);
                if ( cfIndex == (unsigned)-1 ) {
                    Stream*	errorStr = Journal_Register( Error_Type, (Name)self->type  );

                    Journal_Printf( errorStr, "Error- in %s: While parsing "
                                    "definition of wallVC \"%s\" (applies to wall \"%s\"), the cond. func. applied to "
                                    "variable \"%s\" - \"%s\" - wasn't found in the c.f. register.\n",
                                    __func__, self->_dictionaryEntryName, SphereBC_WallEnumToStr[self->_wall],
                                    self->_entryTbl[entry_I].varName, funcName );
                    Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");
                    ConditionFunction_Register_PrintNameOfEachFunc( self->conFunc_Register, errorStr );
                    Journal_Printf( errorStr, ")\n");
                    assert(0);
                }
                self->_entryTbl[entry_I].value.as.typeCFIndex = cfIndex;
            } else if (0 == strcasecmp(valType, "array")) {
                Dictionary_Entry_Value*	valueElement;
                Index			i;

                self->_entryTbl[entry_I].value.type = VC_ValueType_DoubleArray;
                self->_entryTbl[entry_I].value.as.typeArray.size = Dictionary_Entry_Value_GetCount(valueEntry);
                self->_entryTbl[entry_I].value.as.typeArray.array = Memory_Alloc_Array( double,
                        self->_entryTbl[entry_I].value.as.typeArray.size, "SphereBC->_entryTbl[].value.as.typeArray.array" );

                for (i = 0; i < self->_entryTbl[entry_I].value.as.typeArray.size; i++) {
                    valueElement = Dictionary_Entry_Value_GetElement(valueEntry, i);
                    self->_entryTbl[entry_I].value.as.typeArray.array[i] =
                        Dictionary_Entry_Value_AsDouble(valueElement);
                }
            } else if( 0 == strcasecmp( valType, "double" ) || 0 == strcasecmp( valType, "d" ) ||
                       0 == strcasecmp( valType, "float" ) || 0 == strcasecmp( valType, "f" ) ) {
                self->_entryTbl[entry_I].value.type = VC_ValueType_Double;
                self->_entryTbl[entry_I].value.as.typeDouble = Dictionary_Entry_Value_AsDouble( valueEntry );
            } else if( 0 == strcasecmp( valType, "integer" ) || 0 == strcasecmp( valType, "int" ) || 0 == strcasecmp( valType, "i" ) ) {
                self->_entryTbl[entry_I].value.type = VC_ValueType_Int;
                self->_entryTbl[entry_I].value.as.typeInt = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
            } else if( 0 == strcasecmp( valType, "short" ) || 0 == strcasecmp( valType, "s" ) ) {
                self->_entryTbl[entry_I].value.type = VC_ValueType_Short;
                self->_entryTbl[entry_I].value.as.typeShort = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
            } else if( 0 == strcasecmp( valType, "char" ) || 0 == strcasecmp( valType, "c" ) ) {
                self->_entryTbl[entry_I].value.type = VC_ValueType_Char;
                self->_entryTbl[entry_I].value.as.typeChar = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
            } else if( 0 == strcasecmp( valType, "pointer" ) || 0 == strcasecmp( valType, "ptr" ) || 0 == strcasecmp( valType, "p" ) ) {
                self->_entryTbl[entry_I].value.type = VC_ValueType_Ptr;
                self->_entryTbl[entry_I].value.as.typePtr = (void*) ( (ArithPointer)Dictionary_Entry_Value_AsUnsignedInt( valueEntry ));
            } else {
                /* Assume double */
                Journal_DPrintf(
                    Journal_Register( InfoStream_Type, (Name)"myStream"  ),
                    "Type to variable on variable condition not given, assuming double\n" );
                self->_entryTbl[entry_I].value.type = VC_ValueType_Double;
                self->_entryTbl[entry_I].value.as.typeDouble = Dictionary_Entry_Value_AsDouble( valueEntry );
            }
        }
    } else {
        self->_wall = SphereBC_Wall_Size;
        self->_entryCount = 0;
        self->_entryTbl = NULL;
    }
}


void _SphereBC_Delete( void* wallVC )
{
    SphereBC* self = (SphereBC*)wallVC;

    /* Stg_Class_Delete parent */
    _VariableCondition_Delete(self);
}

void _SphereBC_Destroy( void* wallVC, void* data )
{
    SphereBC* self = (SphereBC*)wallVC;

    if (self->_entryTbl) Memory_Free(self->_entryTbl);

    /* Stg_Class_Delete parent */
    _VariableCondition_Destroy( self, data );
}

void* _SphereBC_Copy( void* wallVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap )
{
    SphereBC*		self = (SphereBC*)wallVC;
    SphereBC*		newSphereBC;
    PtrMap*		map = ptrMap;
    Bool		ownMap = False;

    if( !map ) {
        map = PtrMap_New( 10 );
        ownMap = True;
    }

    newSphereBC = (SphereBC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );

    newSphereBC->_dictionaryEntryName = self->_dictionaryEntryName;
    newSphereBC->_wall = self->_wall;
    newSphereBC->_entryCount = self->_entryCount;

    if( deep ) {
        newSphereBC->_mesh = (Mesh*)Stg_Class_Copy( self->_mesh, NULL, deep, nameExt, map );

        if( (newSphereBC->_entryTbl = PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
            newSphereBC->_entryTbl = Memory_Alloc_Array( WallVC_Entry, newSphereBC->_entryCount, "SphereBC->_entryTbl");
            memcpy( newSphereBC->_entryTbl, self->_entryTbl, sizeof(WallVC_Entry) * newSphereBC->_entryCount );
            PtrMap_Append( map, newSphereBC->_entryTbl, self->_entryTbl );
        }
    } else {
        newSphereBC->_mesh = self->_mesh;
        newSphereBC->_entryTbl = self->_entryTbl;
    }

    if( ownMap ) {
        Stg_Class_Delete( map );
    }

    return (void*)newSphereBC;
}


void _SphereBC_Build(  void* wallVC, void* data )
{
    SphereBC*			self = (SphereBC*)wallVC;

    _SphereBC_BuildSelf( self, data );

    _VariableCondition_Build( self, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SphereBC_AssignFromXML( void* wallVC, Stg_ComponentFactory* cf, void* data )
{
}

void _SphereBC_BuildSelf(  void* wallVC, void* data )
{
    SphereBC*			self = (SphereBC*)wallVC;

    if( self->_mesh ) {
        Stg_Component_Build( self->_mesh, data, False );
    }
}


IndexSet* _SphereBC_GetSet(void* variableCondition)
{
    SphereBC* self = (SphereBC*)variableCondition;
    IndexSet* set = NULL;

    switch (self->_wall) {
    case SphereBC_Wall_Inner:
        set = ((SphericalFeMesh*)self->_mesh)->innerSet;
        break;

    case SphereBC_Wall_Outer:
        set = ((SphericalFeMesh*)self->_mesh)->outerSet;
        break;

    case SphereBC_Wall_Size:
    default:
        assert(0);
        break;
    }

    /* Here we make a dynamically allocated copy of the set. 
       The CompositeVC data structure will delete the copy and not the original set
       ( Why the compositeVC wants to delete this data structure i'm not exactly sure - JG Feb14
     */
    return IndexSet_Duplicate(set);
}


VariableCondition_VariableIndex _SphereBC_GetVariableCount(void* variableCondition, Index globalIndex)
{
    SphereBC*	self = (SphereBC*)variableCondition;

    return self->_entryCount;
}


Variable_Index _SphereBC_GetVariableIndex(void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex)
{
    SphereBC*		self = (SphereBC*)variableCondition;
    Variable_Index	searchedIndex = 0;
    Stream*		errorStr = Journal_Register( Error_Type, (Name)self->type  );
    Name		varName;

    varName = self->_entryTbl[varIndex].varName;
    searchedIndex = Variable_Register_GetIndex(self->variable_Register, varName );

    Journal_Firewall( ( searchedIndex < self->variable_Register->count ), errorStr, "Error- in %s: searching for index of "
                      "varIndex %u (\"%s\") at global node number %u failed - register returned index %u, greater than "
                      "count %u.\n", __func__, varIndex, varName, globalIndex, searchedIndex, self->variable_Register->count );

    return searchedIndex;
}


VariableCondition_ValueIndex _SphereBC_GetValueIndex(
    void*				variableCondition,
    Index				globalIndex,
    VariableCondition_VariableIndex	varIndex)
{
    return varIndex;
}


VariableCondition_ValueIndex _SphereBC_GetValueCount(void* variableCondition)
{
    SphereBC*	self = (SphereBC*)variableCondition;

    return self->_entryCount;
}


VariableCondition_Value _SphereBC_GetValue(void* variableCondition, VariableCondition_ValueIndex valIndex)
{
    SphereBC*	self = (SphereBC*)variableCondition;

    return self->_entryTbl[valIndex].value;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/


