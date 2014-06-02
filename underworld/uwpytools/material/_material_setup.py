# Fields package - configures Fields and Variables on a Mesh
import uwpytools as _uw
##############################################################################
## This code adds what is required to the python dictionary 
## to set up Fields and Variables for Underworld.
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary to create materials
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

########################################################################################################################

def materialCreate(componentName="background", rheologyName="", shapeName="", density="1.0"):
    """
    Creates a new material

    Requires a rheology and a shape to create a material with an initial shape.

    Args:
       rheologyName: 
       shapeName:
    """
    globalDict = _uw.GetCurrentPythonDictionary()

    if rheologyName=="":
        _uw.utils.sendError("Must specify a Rheology name")
    if shapeName == "":
        _uw.utils.sendError("Must specify a Shape name")

    _uw.utils.warnMissingComponent(globalDict, shapeName )
    _uw.utils.warnMissingComponent(globalDict, rheologyName )

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)
    
    newComponentDict = _uw.NewComponentEntryInStgDict( globalDict,
                                                       name     = componentName,
                                                       Type     = "RheologyMaterial",
                                                       density  = str(density),
                                                       Shape    = shapeName,
                                                       Rheology = rheologyName
                                                       )

    return newComponentDict
