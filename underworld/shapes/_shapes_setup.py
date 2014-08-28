# Shapes package - configures Shapes
import underworld as _uw
##############################################################################
# This code adds what is required to the python dictionary
# to set up Shapes for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for various Shapes
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''


def everywhereCreate( componentName="backgroundShape"):
    """
    Create the main background "Everywhere" shape.
    """
    globalDict = _uw.dictionary.GetDictionary()

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newComponentShapeDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name     = componentName,
                                                                    Type     = "Everywhere"
                                                                    )

    return newComponentShapeDict


def belowCosinePlaneCreate( componentName="belowCosinePlaneShape",
                            offset="0.2", amplitude="0.02",
                            wavelength="1.8284", phase="0.0"):
    """
    Create the main background "BelowCosinePlane" shape.
    The defaults set up the Rayleigh Taylor benchmark model shape for the lower material
    """
    globalDict = _uw.dictionary.GetDictionary()

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newComponentShapeDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name       = componentName,
                                                                    Type       = "BelowCosinePlane",
                                                                    offset     = str(offset),
                                                                    amplitude  = str(amplitude),
                                                                    wavelength = str(wavelength),
                                                                    phase      = str(phase)
                                                                    )

    return newComponentShapeDict


def boxCreate( componentName="boxShape",
               startList="", endList=""):

    globalDict = _uw.dictionary.GetDictionary()

    sZ = 0.0
    eZ = 1.0
    if startList.__len__() == 3:
        sZ = startList[2]
    if endList.__len__() == 3:
        eZ = endList[2]
    newComponentShapeDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name   = componentName,
                                                                    Type   = "Box",
                                                                    startX = startList[0],
                                                                    startY = startList[1],
                                                                    startZ = sZ,
                                                                    endX   = endList[0],
                                                                    endY   = endList[1],
                                                                    endZ   = eZ
                                                                    )
    return newComponentShapeDict


def sphereCreate( componentName="sphereShape",
                  centreList="", radius="1.0"):

    globalDict = _uw.dictionary.GetDictionary()

    cZ = 0.0
    if centrelist.__len__() == 3:
        cZ = centreList[2]
    newComponentShapeDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name    = componentName,
                                                                    Type    = "Sphere",
                                                                    CentreX = centreList[0],
                                                                    CentreY = centreList[1],
                                                                    CentreZ = cZ,
                                                                    radius  = str(radius)
                                                                    )
    return newComponentShapeDict


def voxelFieldShape(componentName, LowerLimit, UpperLimit, voxel_field="Materials_VoxelField"):
    """
    Creates a voxel field shape from the VoxelField variable

    Args:
        componentName (String): give it a name
        voxel_field (String): name of voxel field in the dictionary (default: "Materials_VoxelField")
        LowerLimit (Float): lower bound
        UpperLimit (Float): upper bound
    Returns:
        None

    """
    globalDict = _uw.dictionary.GetDictionary()
    _uw.utils.warnMissingComponent(globalDict, voxel_field)

    newComponentShapeDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                    name = str(componentName),
                                                                    Type = "FieldValueShape",
                                                                    ValueField = str(voxel_field),
                                                                    LowerLimit = LowerLimit,
                                                                    UpperLimit = UpperLimit
                                                                    )
    return newComponentShapeDict