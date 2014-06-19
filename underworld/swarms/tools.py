# various tools for swarms

from libUnderworld import StGermain

def Swarm_GetVariables( swarm ):
    """
    Returns a python list of swarm variables. 
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        variables (list): python list of variables, variable data stored as tuples: (name, pointer, datatype)
    """
    varList = []
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    for i in range(0,variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes,0)
            # guyTup = (guy.name, guy, SwarmVariable_GetType(datatype))
            guyList = [guy.name, guy, SwarmVariable_GetType(datatype)]
            varList.append(guyList)
            #varList.append(guyTup)
    return varList



def Swarm_GetVariablesAsDict( swarm ):
    """
    Returns a python dictionary of swarm variables with variable names as keys. 
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        variables (dict): python dictionary of variables by variableName, variable data stored as list: (pointer, datatype)
    """

    varDict = {}
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    for i in range(0,variableCount):
        swarmObject = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if swarmObject.variable:
            datatype = StGermain.VariableTypeArrayDeref(swarmObject.variable.dataTypes,0)
            varDict[swarmObject.name] = {}
            varDict[swarmObject.name]['swarmVariable'] = swarmObject
            varDict[swarmObject.name]['dataType'] = SwarmVariable_GetType(datatype)
    
    return varDict

def Swarm_PrintVariables( swarm ):
    """
    Formatted print to standard out of swarm's variables.
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        Nothing
    """
    varNameList = []
    varTypeList = []
    maxLen = 0
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    goodVarCount = 0
    for i in range(0,variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            goodVarCount = goodVarCount + 1
            varNameList.append(str(guy.name))
            if len(guy.name) > maxLen:
                maxLen = len(guy.name)
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes,0)
            varTypeList.append(SwarmVariable_GetType(datatype))
            #print "Name = %-40s Type = %-40s " % (guy.name, SwarmVariable_GetType(datatype))
    for i in range(0,goodVarCount):
        print "Name =", varNameList[i].ljust(maxLen+1), "Type =", varTypeList[i]

def SwarmVariable_GetType( datatype ):
    """
    Returns a string describing the variable type for the swarm variable datatype (as enum) provided
    
    Args:
        datatype (int / enum): swarm variable datatype enumeration
    Returns:
        variableType (str)
    """
    if   datatype == StGermain.Variable_DataType_Variable:
        return "Variable_DataType_Variable"
    elif datatype == StGermain.Variable_DataType_Char:
        return "Variable_DataType_Char"
    elif datatype == StGermain.Variable_DataType_Short:
        return "Variable_DataType_Short"
    elif datatype == StGermain.Variable_DataType_Int:
        return "Variable_DataType_Int"
    elif datatype == StGermain.Variable_DataType_Float:
        return "Variable_DataType_Float"
    elif datatype == StGermain.Variable_DataType_Double:
        return "Variable_DataType_Double"
    elif datatype == StGermain.Variable_DataType_Pointer:
        return "Variable_DataType_Pointer"
    elif datatype == StGermain.Variable_DataType_Size:
        return "Variable_DataType_Size"
    else:
        return "UNKNOWN"


def SwarmVariable_GetValueAt( swarmVar, localParticleIndex ):
    """
    Returns the swarm variable value from a local particle
    
    Args:
        swarmVar (Swig SwarmVariable*): Swarm variable to query.
        localParticleIndex (int)      : Local particle index for particle of interest
    Returns:
        result (list(swarmVar(type))) : The result as a list.  Value(s) in list are of the same type as the swarm variable.
    """
    toreturn = []
    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes,0)
    if   datatype == StGermain.Variable_DataType_Int:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtInt(    swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Float:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtFloat(  swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Double:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtDouble( swarmVar.variable, localParticleIndex, ii ))
    else:
        return "Sorry, the swarm variable datatype", SwarmVariable_GetType(dataType), "is not supported."

    return toreturn


def SwarmVariable_SetValueAt( swarmVar, localParticleIndex, values ):
    """
    Sets the swarm variable value for a local particle
    
    Args:
        swarmVar (Swig SwarmVariable*): Swarm variable to set value for.
        localParticleIndex (int)      : Local particle index for particle of interest
        values (list)                 : values to set on particle
    Returns:
        None
    """
    from libUnderworld import c_arrays
    
    if( swarmVar.dofCount < len(values) ):
        return "Error: size of values list is greater than variable dofCount"

    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes,0)

    if   datatype == StGermain.Variable_DataType_Int:
        valuePtr = c_arrays.IntArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Float:
        valuePtr = c_arrays.FloatArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Double:
        valuePtr = c_arrays.DoubleArray(swarmVar.dofCount)
    else:
        return "Sorry, the swarm variable datatype", SwarmVariable_GetType(dataType), "is not supported."

    for ii in range(0,len(values)):
        valuePtr[ii] = values[ii]

    StGermain.Variable_SetValue( swarmVar.variable, localParticleIndex, valuePtr.cast() )



 
