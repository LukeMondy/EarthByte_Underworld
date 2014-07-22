# various tools for swarms

from libUnderworld import StGermain
import underworld as _uw
import numpy


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
    for i in range(0, variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes, 0)
            # guyTup = (guy.name, guy, SwarmVariable_GetType(datatype))
            guyList = [guy.name, guy, SwarmVariable_GetType(datatype)]
            varList.append(guyList)
            # varList.append(guyTup)
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
    for i in range(0, variableCount):
        swarmObject = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if swarmObject.variable:
            datatype = StGermain.VariableTypeArrayDeref(swarmObject.variable.dataTypes, 0)
            varDict[swarmObject.name] = {}
            varDict[swarmObject.name]['swarmVariable'] = swarmObject
            varDict[swarmObject.name]['dataType'] = SwarmVariable_GetType(datatype)

    return varDict


def Swarm_PrintVariables( swarmIn ):
    """
    Formatted print to standard out of swarm's variables.
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.

    Args:
        swarm (Swig Swarm*, str): The swarm, either provided as a swig generated pointer, or via the textual name.
    Returns:
        Nothing
    """
    if isinstance(swarmIn,(str)):
        swarm = _uw._stgermain.GetLiveComponent(swarmIn)
    else:
        swarm = swarmIn

    varNameList = []
    varTypeList = []
    maxLen = 0
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    goodVarCount = 0
    for i in range(0, variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            goodVarCount = goodVarCount + 1
            varNameList.append(str(guy.name))
            if len(guy.name) > maxLen:
                maxLen = len(guy.name)
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes, 0)
            varTypeList.append(SwarmVariable_GetType(datatype))
            # print "Name = %-40s Type = %-40s " % (guy.name, SwarmVariable_GetType(datatype))
    for i in range(0, goodVarCount):
        print "Name =", varNameList[i].ljust(maxLen + 1), "Type =", varTypeList[i]


def SwarmVariable_GetType( datatype ):
    """
    Returns a string describing the variable type for the swarm variable datatype (as enum) provided

    Args:
        datatype (int / enum): swarm variable datatype enumeration
    Returns:
        variableType (str)
    """
    if datatype == StGermain.Variable_DataType_Variable:
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
    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes, 0)
    if datatype == StGermain.Variable_DataType_Int:
        for ii in range(0, swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtInt(    swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Float:
        for ii in range(0, swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtFloat(  swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Double:
        for ii in range(0, swarmVar.dofCount):
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

    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes, 0)

    if datatype == StGermain.Variable_DataType_Int:
        valuePtr = c_arrays.IntArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Float:
        valuePtr = c_arrays.FloatArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Double:
        valuePtr = c_arrays.DoubleArray(swarmVar.dofCount)
    else:
        return "Sorry, the swarm variable datatype", SwarmVariable_GetType(dataType), "is not supported."

    for ii in range(0, len(values)):
        valuePtr[ii] = values[ii]

    StGermain.Variable_SetValue( swarmVar.variable, localParticleIndex, valuePtr.cast() )


def ArrayParticleLayoutCreate(  arrParticles, componentName = "manualLayout", dim=2):
  """
  This function sets up a manual particle layout using the points from an array. 
  This works best if arrParticles is a numpy array / matrix of nParticles x dim dimensions.

  The 'ManualParticleLayout' component python object is returned.

"""
  nParticles = len(arrParticles)
  manualPoints = [dict() for x in range(nParticles)]


  for m in range(0,nParticles):
    point = arrParticles[m]
    manualPoints[m]["x"] = str(point[0])
    manualPoints[m]["y"] = str(point[1])
    if dim == 3:
      manualPoints[m]["z"] = str(point[2])
    

  newLayout = _uw.dictionary.UpdateDictWithComponent( name   = componentName,
                                              Type   = "ManualParticleLayout",
                                              manualParticlePositions = manualPoints
                                              )
  

  return newLayout

def TracerSwarmCreate( particleLayout, 
                      componentName="tracerSwarm", 
                      meshName=""
                      ):
    """
    Create a swarm intended for use as tracers
    """

    tracerSwarm = _uw.swarms.setup.materialSwarmCreate(componentName=componentName,
                                meshName = meshName,
                                particleLayout=particleLayout)
    _uw.swarms.setup.materialSwarmAdvectorCreate(componentName=componentName + "-Advector",
                                        swarmName=componentName,
                                        velocityFieldName = "VelocityField"
                                        )

   

    return tracerSwarm

def InterpolateSwarmVariable(swarm,fieldname,dim=2):

    """
      This function interpolates the given field over each point position in a swarm. 
      It returns a dictionary, with a 'position' array and a 'field' array.
      
      Args:
        
        swarm:       a python swarm object
                            

        fieldname:   string -  the name of the field that will be interpolated
                             

        dim:         model dimension, defaults to 2


    """

    
    tracerSwarm = _uw._stgermain.GetLiveComponent(swarm)
    tracerVariables = Swarm_GetVariablesAsDict(tracerSwarm)
    nParticles = tracerSwarm.particleLocalCount
    
    # Get positions of swarm particles

    pos = numpy.zeros((nParticles,dim))

    for j in range(0,nParticles):
        pos[j][0]=SwarmVariable_GetValueAt(tracerVariables[swarm + "-PositionX"]["swarmVariable"],j)[0]
        pos[j][1]=SwarmVariable_GetValueAt(tracerVariables[swarm + "-PositionY"]["swarmVariable"],j)[0]
        if dim==3:
            pos[j][2]=SwarmVariable_GetValueAt(tracerVariables[swarm + "-PositionZ"]["swarmVariable"],j)[0]
    
    # Interpolate field at each position

    arrField = numpy.zeros((len(pos),dim))  
    
    for j in range(0,len(pos)):
        fieldVal = _uw.fields.tools.FieldVariable_InterpolateValueAt( fieldname, tuple(pos[j]) )  
        arrField[j][0] = list(fieldVal)[0][0]
        arrField[j][1] = list(fieldVal)[0][1]
        if dim==3:
            arrField[j][2] = list(fieldVal)[0][2]
    

    return {'position' : pos, 'field' : arrField}

