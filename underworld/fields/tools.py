# Various tools for interrogating Underworld fields

from libUnderworld import StgDomain
from libUnderworld import StgFEM
import underworld._stgermain as _stgermain

def FieldVariable_InterpolateValueAt( fieldVar, coord ):
    """
    Returns the interpolated result of a field variable.

    Args:
        fieldVar (str, Swig FieldVariable*): Field variable to query, either provided by name (str) or with a swig pointer.
        coord (tuple(float)): The location to query
    Returns:
        result (list(float)): The interpolated result as a list
        success (str): "OTHER_PROC"     - Query location is on another processor
                       "LOCAL"          - Query location is on this processor
                       "SHADOW"         - Query location is in shadow space of this processor
                       "OUTSIDE_GLOBAL" - Query location is outside the global domain of this field variable
    """
    from libUnderworld import c_arrays

    if type(fieldVar)==str:
        fieldVar = _stgermain.GetLiveComponent(fieldVar)

    result = c_arrays.DoubleArray(fieldVar.fieldComponentCount)

    InterpolationResult = StgDomain.FieldVariable_InterpolateValueAt( fieldVar, coord, result.cast() )

    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"

    toreturn = []
    for i in range(0, fieldVar.fieldComponentCount):
        toreturn.append(result.__getitem__(i))

    return toreturn, InterpolationResult


def FieldVariable_GetMinFieldMagnitude( fieldVar ):
    """
    Returns the minimum value of a field variable

    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (float): the minimum field variable magnitude
    """
    return StgDomain.FieldVariable_GetMinGlobalFieldMagnitude( fieldVar )


def FieldVariable_GetMaxFieldMagnitude( fieldVar ):
    """
    Returns the maximum value of a field variable

    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        max (float): the maximum field variable magnitude
    """
    return StgDomain.FieldVariable_GetMaxGlobalFieldMagnitude( fieldVar )


def FieldVariable_GetMinAndMaxLocalCoords( fieldVar ):
    """
    Returns the domain of this field variable local to this process as two tuples (min,max)

    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (tuple(float)): the minimum local coordinate range
        max (tuple(float)): the maximum local coordinate range
    """

    result = StgDomain.FieldVariable_GetMinAndMaxLocalCoords( fieldVar )
    return result[0], result[1]


def FieldVariable_GetMinAndMaxGlobalCoords( fieldVar ):
    """
    Returns the domain of this field variable globally as two tuples (min,max)

    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (tuple(float)): the minimum global coordinate range
        max (tuple(float)): the maximum global coordinate range
    """

    result = StgDomain.FieldVariable_GetMinAndMaxGlobalCoords( fieldVar )
    return result[0], result[1]


def FeVariable_Integrate( feVar, gaussSwarm=None ):
    """
    Returns the integral of this FE variable.
    Note that this is currently only compatible with scalar fields.

    Args:
        fieldVar (Swig FeVariable*): Field Variable to query
        gaussSwarm (Swig Swarm*): (optional) Gauss swarm to integrate over. If nothing provided,
                                             will try and locate gauss swarm with name "gaussSwarm"
    Returns:
        integral (float): the integral of the FeVariable.
    """
    if gaussSwarm == None:
        gaussSwarm = _stgermain.GetLiveComponent("gaussSwarm")

    if feVar.fieldComponentCount != 1:
        print "Error: This function currently only supports scalar fields."
        return

    return StgFEM.FeVariable_Integrate( feVar, gaussSwarm )
