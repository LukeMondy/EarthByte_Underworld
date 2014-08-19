# Rheology package
import underworld as _uw
##############################################################################
# This code adds what is required to the python dictionary
# to set up various Rheologies for Underworld.
# We eventually pass the python dictionary back to Underworld
# and Underworld then uses this information to configure and set
# itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for various Rheologies
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''


def isoviscousCreate( componentName="backgroundViscosity", eta0="eta0"):
    """
    Create an isoviscous constant viscosity rheology.
    """
    globalDict = _uw.dictionary.GetDictionary()

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newComponentIsoviscosityDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                           name = componentName,
                                                                           Type = "MaterialViscosity",
                                                                           eta0 = str(eta0)
                                                                           )

    return newComponentIsoviscosityDict


def vonMisesCreate(componentName="vonMisesYieldRheology", eta0="1.0"):
    """
    Set up a Von Mises Yield rheology.
    """

    globalDict = _uw.dictionary.GetDictionary()

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                      name = componentName,
                                                      Type = "VonMises",
                                                      eta0 = str(eta0)
                                                      )
    return newDict


def arrheniusCreate(componentName="arrheniusRheology", eta0="1.0e-6", temperatureField="TemperatureField", activationEnergy="27.63102112"):
    """
    Set up  Arrhenius (temperature dependent) rheology
    """

    globalDict = _uw.dictionary.GetDictionary()

    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    newDict = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                      name = componentName,
                                                      Type = "Arrhenius",
                                                      eta0 = str(eta0),
                                                      TemperatureField  = temperatureField,
                                                      activationEnergy = str(activationEnergy)
                                                      )
    return newDict


# should this one live in geometry?
def joinRheologyAndShape(componentName="background", rheologyName="", shapeName="", density="1.0"):
    """
    Set the Rheology associated with a Shape and vice-versa
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
                                                               name     = componentName,
                                                               Type     = "RheologyMaterial",
                                                               density  = str(density),
                                                               Shape    = shapeName,
                                                               Rheology = rheologyName
                                                               )

    return newComponentDict


def TemperatureDependentConductivity(thermalConductivity, componentName="", a=0.5):
  """
  Temperature-dependent thermal conductivity (k)
  k = k0 * (298/T)^a     where k0 is room temp conductivity, and 0 < a < 1

  Similar to arrheniusRheology but for geothermal use

  Args:
    componentName (String): give it a name
    thermalConductivity (Float): k0
    a (Float): power function, a
  """
  globalDict = _uw.dictionary.GetDictionary()

  # Make sure temperature field exists
  if "temperatureField" not in globalDict | "TemperatureField" not in globalDict:
    temperatureField = globalDict["info"]["TemperatureField"]
  # Make sure the solver is non-linear
  if "energySLE" not in globalDict:
    energySolver = _uw.dictionary.UpdateDictWithComponent(  globalDict,
                                                            Type = "Energy_SLE_Solver"
                                                            )
    energySLE = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                        name = "energySLE",
                                                        Type = "Energy_SLE",
                                                        SLE_Solver = "energySolver",
                                                        Context = "context",
                                                        StiffnessMatrix = "kMatrix",
                                                        ForceVector = "fVector",
                                                        SolutionVector = "sVector",
                                                        isNonLinear = "yes"
                                                        )
  if "temperatureFieldPpc" not in globalDict:
    temperatureFieldPpc = _uw.dictionary.UpdateDictWithComponent( globalDict,
                                                                  Type = "Ppc_Variable",
                                                                  FieldVariable = "TemperatureField"
                                                                  )

  #componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

  newComponentDict = uw.dictionary.UpdateDictWithComponent( globalDict,
                                                            name = componentName,
                                                            Type = "Ppc_Polynomial",
                                                            ImposeMinMax = True,
                                                            MaxValue = 5.0,
                                                            MinValue = 0.1,
                                                            Field = "TemperatureField",
                                                            Terms = {"Coefficient" : float(thermalConductivity)*298**a, "Power" : -a}
                                                            )
  return newComponentDict