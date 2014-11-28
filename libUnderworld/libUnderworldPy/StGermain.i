/* -*- C -*-  (not really, but good for syntax highlighting) */

%module StGermain

%{
/* Includes the header in the wrapper code */
#include <StGermain/StGermain.h>
%}


%include "StGermain_Typemaps.i"

int    Variable_GetValueAtInt(Variable*   , unsigned int, unsigned int);
float  Variable_GetValueAtFloat(Variable* , unsigned int, unsigned int);
double Variable_GetValueAtDouble(Variable*, unsigned int, unsigned int);

%include "libStGermain/main.h"
%include "Base/Automation/types.h"
%include "Base/Context/types.h"
%include "Base/Extensibility/types.h"
%include "Base/Foundation/types.h"
%include "Base/IO/types.h"
%include "Utils/types.h"
%include "Base/Foundation/Class.h"
%include "Base/Foundation/Object.h"
%include "Base/Automation/Stg_Component.h"
%include "Base/Automation/LiveComponentRegister.h"       
%include "Base/Automation/Stg_ComponentFactory.h"       
%include "Base/Automation/Stg_ComponentRegister.h"
%include "Base/Context/AnalyticFunction.h"
%include "Base/Context/Variable.h"
%include "Base/Context/AnalyticFunction_Register.h"
%include "Base/Context/ConditionFunction.h"       
%include "Base/Context/ConditionFunction_Register.h"
%include "Base/Context/Variable_Register.h"
%include "Base/Context/VariableCondition_Register.h"
%include "Base/Context/VariableCondition.h"
%include "Base/Extensibility/EntryPoint.h"
%include "Base/Extensibility/Init.h"
%include "Base/Extensibility/EntryPoint_Register.h"      
%include "Base/Extensibility/ExtensionManager_Register.h"      
%include "Base/Extensibility/ModulesManager.h"
%include "Base/Foundation/NamedObject_Register.h"       
%include "Base/Container/HashTable.h"
%include "Base/Automation/HierarchyTable.h"
%include "Base/IO/IO_Handler.h"       
%include "Base/IO/XML_IO_Handler.h"       
%include "Base/IO/File.h"
%include "Base/IO/Dictionary.h"
%include "Base/IO/DictionaryUtils.h"
%include "Base/IO/Dictionary_Entry.h"
%include "Base/IO/Dictionary_Entry_Value.h"
%include "Base/Foundation/ObjectList.h"       
%include "Base/Foundation/ObjectList.h"
%include "Base/Context/AbstractContext.h"
%include "Base/Context/Codelet.h"
%include "Utils/DummyComponent.h"
%include "Utils/Scaling.h"
