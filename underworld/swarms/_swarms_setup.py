# Swarms package - configures swarms
import underworld as _uw


##############################################################################
## This code adds what is required to the python dictionary 
## to set up a Swarm for Underworld.
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for an Integration Swarm. It can be either particles or Gauss points.
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

#    A PIC integration Swarm needs a Mapper.
#    The Mapper in turn needs an integration Swarm and a material Swarm.
#    So best to create Mapper and integration Swarm simultaneously.
#    This function requires an existing material Swarm to exist in the Dictionary.

def integrationSwarmCreate(  componentName="",
                             swarmType   = "PIC",
                             meshName    = "",              # Maybe should warn if Mesh actually doesn't exist yet?
                             integratorName    = "timeIntegrator",
                             weightsName       = "weights",
                             materialSwarmName = "",
                             mapperType        = "CoincidentMapper",
                             contextName       = "context",             # Maybe should warn if context actually doesn't exist yet?
                             auto              = False
                             ):

    """
    Set up an integration Swarm in the dictionary.
    Possible swarmTypes are ( PIC, Gauss )
    
    When swarmType is PIC:
       An integration Swarm needs a Mapper.
       The Mapper in turn needs the integration Swarm and a material Swarm.
       So this function requires an existing material Swarm to exist in the Dictionary
       and we will create a Mapper at same time as this integration Swarm.


    Args:
        componentName (string)         : The name of the Swarm (optional: will create a name based on swarmType if null)
        swarmType (string)             : Must be one of ( "PIC", "pic", "Gauss", "gauss" ). Defaults to PIC
        meshName (string)              : This swarm needs a meshName (this is overridden if swarmType is PIC)
        integratorName (string)        : Will create this if does not exist already.
        particlesPerCell int           : If using PIC then set number of initial particles per cell
        weightsName (string)           : Required for PIC. Will create this if does not exist already.
        materialSwarmName (string)     : Required for PIC. Must already exist
        mapperType                     : Required for PIC.

    Returns:

        
    """

    swarmAllowedTypes=["PIC", "pic", "Gauss", "gauss"]
    
    # mcount=0
    # for m in swarmAllowedTypes:
    #   if( swarmType != m ):
    #         mcount += 1

    # if (not mcount):
    #   print "Error: Need to choose one of  ( PIC, Gauss ) for swarmType"
    #   return
    if meshName=="":
        print "Error: Need to provide a meshName"
        return
    
    if not swarmType in swarmAllowedTypes:
      print "Error: Need to choose one of  ( PIC, Gauss ) for swarmType"
      return

    globalDict = _uw.dictionary.GetDictionary()

    if swarmType in ["PIC", "pic"]:
        # we have a chicken and egg scenario here.
        # To create a mapper need an integration point swarm and vice-versa
        # So we will create both at same time here
        # But then we need a material Swarm.

        if materialSwarmName == "" or materialSwarmName not in globalDict["components"]:
            print "Error: Need a material Swarm to exist and it's name passed in"
            print "       e.g. materialSwarmName=\"materialSwarm\""
            return -1

        #Now get the Mesh from the material Swarm
        materialSwarm = globalDict["components"][materialSwarmName]
        meshName = materialSwarm["FeMesh"]
        intMapperName = "integrationPointsMapper"
        
    if componentName == "":
        componentName = "integration"+swarmType+"Swarm"    # e.g. integrationPICSwarm
    
    # Need to test if we already have a Swarm of the name we are going to use here.
    componentName = _uw.utils.checkForNewComponentName(globalDict, componentName)

    # At this point we should have ensured a unique name for our new Swarm component
    
    ############################################################################################
    # We need a time Integrator for the swarm.
    # The usual name is "timeIntegrator" so will check for that
    # If it doesn't already exist then let's make one
    ############################################################################################    

    if integratorName not in globalDict["components"]:
        integrator = timeIntegratorCreate(componentName="timeIntegrator", contextName=contextName, order="2", simultaneousFlag="False")
        integratorName = integrator["name"]

    ############################################################################################
    # We need a few other components as well.
    # Just stomping on stuff like Godzilla here...
    ############################################################################################
    ########################################################################################################################################
    # GAUSS integration swarm
    ########################################################################################################################################
 
    if swarmType in ["Gauss", "gauss"]:
        cellLayout     = _uw.dictionary.UpdateDictWithComponent( name="cellLayout", 
                                                       Type="SingleCellLayout")
        particleLayout = _uw.dictionary.UpdateDictWithComponent( name="gaussParticleLayout",
                                                       Type="GaussParticleLayout",
                                                       gaussParticlesX = "gaussParticlesX",
                                                       gaussParticlesY = "gaussParticlesY",
                                                       gaussParticlesZ = "gaussParticlesZ",
                                                      )

        newComponentSwarmDict = _uw.dictionary.UpdateDictWithComponent( name           =            componentName,
                                                              Type           = "IntegrationPointsSwarm",
                                                              CellLayout     =       cellLayout["name"],
                                                              ParticleLayout =   particleLayout["name"],
                                                              FeMesh         =                 meshName,
                                                              TimeIntegrator =           integratorName
                                                                )
    ########################################################################################################################################
    # PIC integration Swarm
    ########################################################################################################################################

    mapper={}

    if swarmType in ["PIC", "pic"]:
        # Now we need some pic specific stuff
        particleLayout    = _uw.dictionary.UpdateDictWithComponent( name="mappedParticleLayout", 
                                                          Type="MappedParticleLayout")
        #elementCellLayout = _uw.dictionary.UpdateDictWithComponent( name="elementCellLayout", Type="ElementCellLayout", Mesh=meshName)
        elementCellLayoutName = materialSwarm["CellLayout"]

        if weightsName not in globalDict["components"]:
            weights = weightCalculatorCreate( componentName=weightsName,
                                               materialPointsSwarmName=materialSwarmName) 
                                               # create a weights component with basic defaults.
            weightsName = weights["name"]

        newComponentSwarmDict = _uw.dictionary.UpdateDictWithComponent( name           =             componentName,
                                                              Type           =  "IntegrationPointsSwarm",
                                                              CellLayout     =     elementCellLayoutName,
                                                              ParticleLayout =    particleLayout["name"],
                                                              FeMesh         =                  meshName,
                                                              TimeIntegrator         =    integratorName,
                                                              IntegrationPointMapper =     intMapperName,
                                                              WeightsCalculator      =       weightsName
                                                            )
 
        # There is some circularity here. We might want different mappers but will do this for now.
        # int swarm needs a mapper; mapper needs an int swarm (as well as a material swarm)
        # might be better to pass in existing mapperName and make a new one if not exist as with the timeIntegrator
        # The mapper also needs a material point swarm.

        matPointSwarmName = newComponentSwarmDict["name"]

        mapper = _uw.dictionary.UpdateDictWithComponent(   name                   =                 intMapperName,
                                                 Type                   =                    mapperType, 
                                                 IntegrationPointsSwarm = newComponentSwarmDict["name"], 
                                                 MaterialPointsSwarm    =             materialSwarmName)
    ########################################################################################################################################
 
    return [newComponentSwarmDict,integratorName, mapper]

############################################################################################################################################

def weightCalculatorCreate( componentName="weights",
                            materialPointsSwarmName="materialSwarm",
                            resX=10, resY=10, resZ=10,
                            lowerT=0.6,
                            upperT=25,
                            maxDeletions=3,
                            maxSplits=3,
                            InFlow=False):
    """
    Create the weights Calculator (This one is used in PIC)
    """

    weights = _uw.dictionary.UpdateDictWithComponent( name="weights",
                                            Type="PCDVC",
                                            MaterialPointsSwarm=materialPointsSwarmName,
                                            resolutionX=resX,
                                            resolutionY=resY,
                                            resolutionZ=resZ,
                                            lowerT=lowerT,
                                            upperT=upperT,
                                            maxDeletions=maxDeletions,
                                            maxSplits=maxSplits,
                                            InFlow=InFlow
                                          )

    return weights


def materialSwarmCreate( componentName="materialSwarm", 
                                           meshName="", 
                                           particlesPerCell="particlesPerCell",
                                           pic=True):
    """
    Create a material points Swarm
    """

    globalDict = _uw.dictionary.GetDictionary()

    if meshName == "" or meshName not in globalDict["components"]:
        print "Error: Need a Mesh to exist and it's name passed in"
        print "       e.g. meshName=\"velocityMesh\""
        return -1

    if pic:
        particleLayout    = _uw.dictionary.UpdateDictWithComponent( name="spaceFillerParticleLayout", 
                                                          Type="SpaceFillerParticleLayout",
                                                          averageInitialParticlesPerCell=particlesPerCell )
    else:
        particleLayout    = _uw.dictionary.UpdateDictWithComponent( name="backgroundParticleLayout", 
                                                                                             Type="BackgroundParticleLayout" )

    elementCellLayout = _uw.dictionary.UpdateDictWithComponent( name="elementCellLayout", 
                                                      Type="ElementCellLayout", 
                                                      Mesh=meshName
                                                    )
    pMovementHandler  = _uw.dictionary.UpdateDictWithComponent( name="pMovementHandler", 
                                                      Type="ParticleMovementHandler"
                                                    )
    pShadowSync       = _uw.dictionary.UpdateDictWithComponent( name="pShadowSync", 
                                                      Type="ParticleShadowSync"
                                                    )

    pList  = [pMovementHandler["name"],pShadowSync["name"]]

    newComponentMaterialSwarmDict = _uw.dictionary.UpdateDictWithComponent( name             =             componentName,
                                                                  Type             =     "MaterialPointsSwarm",
                                                                  CellLayout       = elementCellLayout["name"],
                                                                  ParticleLayout   =    particleLayout["name"],
                                                                  FeMesh           =                  meshName,
                                                                  ParticleCommHandlers = pList,
                                                                  SplittingRoutine = "splittingRoutine",  # this one not needed anymore?
                                                                  RemovalRoutine   = "removalRoutine", 
                                                                  EscapedRoutine   = "escapedRoutine"
                                                                )
    
    return newComponentMaterialSwarmDict


def timeIntegratorCreate( componentName="timeIntegrator", 
                          contextName="context", 
                          order="2", 
                          simultaneousFlag="False"):
    """
    Create a time integrator (needed by integration swarm)
    """
    # todo check for existence of context or just warn if not exist...
    
    integrator = _uw.dictionary.UpdateDictWithComponent( name  = componentName,
                                               Type  = "TimeIntegrator",
                                               order = order,
                                               simultaneous = simultaneousFlag,
                                               Context = contextName)
    return integrator


def materialSwarmAdvectorCreate(componentName  = "materialSwarmAdvector",
                                integratorName = "timeIntegrator",
                                swarmName      = "materialSwarm",
                                velocityFieldName = "velocityField",
                                periodicBCsManagerName = "periodicBCsManager",  # by default no periodic BCs
                                fallBackToFirstOrder   = "True"
                                ):
    """
    Used for PIC
    Create a material Swarm Advector Component in the global Dictionary
    """
    
    newComponentMaterialSwarmAdvectorDict = _uw.dictionary.UpdateDictWithComponent( name  = componentName,
                                                                          Type  = "SwarmAdvector",
                                                                          Swarm = swarmName,
                                                                          TimeIntegrator = integratorName,
                                                                          VelocityField  =velocityFieldName,
                                                                          PeriodicBCsManager = periodicBCsManagerName,
                                                                          allowFallbackToFirstOrder = fallBackToFirstOrder
                                                                          )
    return newComponentMaterialSwarmAdvectorDict
