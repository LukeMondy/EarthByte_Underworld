# Fields package - configures Fields and Variables on a Mesh
import underworld as _uw
##############################################################################
# This code adds what is required to the python dictionary
# to set up Fields and Variables for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary to create materials
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

########################################################################################################################


# maybe let shapeName be a list as well.
# if so then make a union shape from the list of shapes.
# use type(List) is list
# type("jfah") is str
# or can make a union first and just pass in this name
def materialCreate(componentName="background", rheologyName="", shapeName="", density="1.0", alpha="0.0", referenceTemperature="0.0"):
    """
    Creates a new material

    Requires a rheology and a shape to create a material with an initial shape.

    Args:
       rheologyName:
       shapeName:
    """
    globalDict = _uw.dictionary.GetDictionary()

    if rheologyName == "":
        _uw.utils.sendError("Must specify a Rheology name")
    if shapeName == "":
        _uw.utils.sendError("Must specify a Shape name")

    _uw.utils.warnMissingComponent(globalDict, shapeName )
    _uw.utils.warnMissingComponent(globalDict, rheologyName )

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newComponentDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                               name                 = componentName,
                                                               Type                 = "RheologyMaterial",
                                                               density              = str(density),
                                                               Shape                = shapeName,
                                                               Rheology             = rheologyName,
                                                               alpha                = str(alpha),
                                                               referenceTemperature = str(referenceTemperature)

                                                               )

    return newComponentDict


def geothermalMaterialCreate(componentName="", shapeName="", thermalConductivity="", heatProduction="", Tdependence = False, a = 0.5):
  """
  Creates a new material with a geothermal rheology

  Requires a shape to already be defined
  If the rheology should be temperature dependent, set Tdependence to True

  Args:
    componentName (String): give it a name
    shapeName (String): Point to a predefined shape
    thermalConductivity (Float)
    heatProduction (Float)
  """
  globalDict = _uw.dictionary.GetDictionary()

  if componentName == "":
    _uw.utils.sendError("Must specify a component name")
  if shapeName == "":
    _uw.utils.sendError("Must specify a Shape name")
  if thermalConductivity == "":
    _uw.utils.sendError("Must specify a thermal conductivity value")
  if heatProduction == "":
    _uw.utils.sendError("Must specify a heat production value")

  _uw.utils.warnMissingComponent(globalDict, shapeName )
  componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

  if Tdependence is False:
    # thermal conductivity IS NOT temperature dependent
    newComponentDict = _uw.dictionary.UpdateDictWithComponent(  globalDict,
                                                                name = componentName,
                                                                Type = "Material",
                                                                Shape = shapeName,
                                                                thermalConductivity = thermalConductivity,
                                                                heatProduction = heatProduction
                                                                )
    return newComponentDict

  elif Tdependence is True:
    # thermal conductivity IS temperature dependent
    # k = k0 * (298/T)^a     where a is between 0 and 1
    materialPpc = uw.rheology.TemperatureDependentConductivity(str(componentName)+"_MaterialConductivityPpc", float(thermalConductivity), a)

    newComponentDict = uw.dictionary.UpdateDictWithComponent( globalDict,
                                                              name = componentName,
                                                              Type = "Material",
                                                              Shape = shapeName,
                                                              thermalConductivity = str(componentName)+"_MaterialConductivityPpc",
                                                              heatProduction = heatProduction
                                                              )
    return newComponentDict, materialPpc