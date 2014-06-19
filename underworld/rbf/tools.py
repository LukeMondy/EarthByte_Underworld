# Various tools for RBF

# This is the "proper" way to do things ... differentiate the kernel etc etc.
# But this is not particularly stable and needs a lot of points to work.
# We probably need a more sophisticated weighting scheme to make this work
# But calculating a FD gradient is quite accurate ... 

def RBFGradientVariable_InterpolateValueAt1( RBFfieldVar, coord , axis ):
    """
    Returns the interpolated result of an RBF variable.
    
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
    
    import c_arrays
    
    if type(RBFfieldVar)==str:
        RBFfieldVar = GetLiveComponent(RBFfieldVar)

    result = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    coordArray = c_arrays.DoubleArray(3)

    coordArray[0] = coord[0]
    coordArray[1] = coord[1]
    coordArray[2] = coord[2]

    InterpolationResult = Underworld.RBFFieldVariable_InterpolateGradientValueAt( RBFfieldVar, coordArray.cast(), result.cast(), axis )

    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"

    toreturn = []
    for i in range(0,RBFfieldVar.fieldComponentCount):
        toreturn.append(result.__getitem__(i))

    return toreturn, InterpolationResult


# 
# This is a Finite Difference calculation of the gradient of the RBF field 

def RBFGradientVariable_InterpolateValueAt( RBFfieldVar, coord , axis ):
    """
    Returns the interpolated result of an RBF variable.
    
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
    
    import c_arrays
    
    if type(RBFfieldVar)==str:
        RBFfieldVar = GetLiveComponent(RBFfieldVar)

    result0 = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    resultM = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    resultP = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)

    coordArray = c_arrays.DoubleArray(3)
    coordArray[0] = coord[0]
    coordArray[1] = coord[1]
    # if()
    coordArray[2] = coord[2]

    # Compare it with this result

    supportSize = RBFfieldVar.rbfManager.particleSupportRadius
    
    InterpolationResult = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), result0.cast())

    if InterpolationResult == 3:
        return [], "OUTSIDE_GLOBAL"

    coordArray[axis] -= supportSize * 0.5
    InterpolationResultM = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), resultM.cast())

    coordArray[axis] += supportSize
    InterpolationResultP = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), resultP.cast())

    
    # Weights 

    WM = -1.0; W0 =  0.0; WP =  1.0
    delta = supportSize

    if InterpolationResultM == 3:
        WM = 0.0
        W0 = -1.0
        delta = 0.5 * supportSize

    if InterpolationResultP == 3:
        WP = 0.0
        W0 += 1.0  # If both the sample points are outside the domain, the result is zero ... 
        delta = 0.5 * supportSize

    # if(axis == 1):
    #     print "Coord P - {} - {}".format(coordArray[1],coord[1])
    #     print "Weights {}/{}/{} Del {}".format(WM,W0,WP,delta)
    #     print "Values  {}/{}/{} ".format(resultM[0],result0[0],resultP[0])
    #     print "Results {}/{}/{} ".format(InterpolationResultM,InterpolationResult,InterpolationResultP)


    toreturn = []
    for i in range(0,RBFfieldVar.fieldComponentCount):
        toreturn.append( (WM * resultM[i] + W0 * result0[i] + WP * resultP[i]) / delta )


    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"


    return toreturn, InterpolationResult


