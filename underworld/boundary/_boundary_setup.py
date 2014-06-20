# Boundary package - activates and configures boundary conditions
import underworld as _uw
##############################################################################
## This code adds what is required to the python dictionary
## to set up a Boundary Conditions for Underworld.
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for Boundary Conditions
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

def emptyBoundaryConditionCreate(bcEntry=""):
    """
    Create a new empty BC/IC entry in the Dictionary.
    Note that it is not in the components sub-dictionary
    but has it's own entry at the top level.
    (usually "velocityBCs" or  "temperatureBCs" )
    """
    globalDict = _uw.dictionary.GetDictionary()

    if bcEntry=="":
        field = globalDict["info"]["velocityField"]
        bcEntry = globalDict["info"][field+"BCs"]

    newList=[]
    globalDict[bcEntry]=dict()
    globalDict[bcEntry]["name"]=bcEntry
    globalDict[bcEntry]["Type"]="CompositeVC"
    globalDict[bcEntry]["vcList"]=[]

    return

def wallFreeSlipCreate(bcEntry="", wall="left"):
    """
    Set a wall to be Free Slip for a Cartesian Mesh
    """

    globalDict = _uw.dictionary.GetDictionary()

    # Is it a good idea to make sure we have things like this on the dictionary always?
    # it makes some things a lot easier though
    if bcEntry=="":
        field = globalDict["info"]["velocityField"]
        bcEntry = globalDict["info"][field+"BCs"]

    if "dim" in globalDict.keys():
        dim = globalDict["dim"]
        #_uw.utils.sendWarning("Dim set to"+" "+str(dim))
    else:
        dim = 2    # If it's not set I am going to assume it's 2D
        _uw.utils.sendWarning("The value 'dim' not set in dictionary. So assuming dim=2")

    if bcEntry not in globalDict:
        # quiet for now
        #print "No boundary condition entry found with that name"
        #print "So creating one!"
        newBCEntry = emptyBoundaryConditionCreate(bcEntry=bcEntry)
    if "vcList" not in globalDict[bcEntry]:
        globalDict[bcEntry]["vcList"]=[]  # quietly create the necessary list

    # structure here is dict->list->dict->list
    vcListList = globalDict[bcEntry]["vcList"] # get list from BC dictionary (which is attached to main dictionary)

    wallVCdict = dict()
    wallVCdict["type"]="WallVC"
    wallVCdict["wall"]=wall
    wallVCdict["variables"]=[]


    vx=dict()
    vx["name"]="vx"
    vx["type"]="double"
    vx["value"]=0.0
    vy=dict()
    vy["name"]="vy"
    vy["type"]="double"
    vy["value"]=0.0
    vz=dict()
    vz["name"]="vz"
    vz["type"]="double"
    vz["value"]=0.0

    if wall=="left":
        wallVCdict["variables"].append(vx)
    if wall=="right":
        wallVCdict["variables"].append(vx)
    if wall=="top":
        wallVCdict["variables"].append(vy)
    if wall=="bottom":
        wallVCdict["variables"].append(vy)

    if dim == 3:
        if wall=="front":
            wallVCdict["variables"].append(vz)
        if wall=="back":
            wallVCdict["variables"].append(vz)

    vcListList.append(wallVCdict)

    return globalDict[bcEntry]


def wallNoSlipCreate(bcEntry="", wall="left"):
    """
    Set a wall to be No Slip for a Cartesian Mesh
    """

    globalDict = _uw.dictionary.GetDictionary()

    # Is it a good idea to make sure we have things like this on the dictionary always?
    # it makes some things a lot easier though
    if bcEntry=="":
        field = globalDict["info"]["velocityField"]
        bcEntry = globalDict["info"][field+"BCs"]

    if "dim" in globalDict.keys():
        dim = globalDict["dim"]
    else:
        dim = 2    # If it's not set I am going to assume it's 2D
        _uw.utils.sendWarning("The value 'dim' not set in dictionary. So assuming dim=2")

    if bcEntry not in globalDict:
        #print "No boundary condition entry found with that name"
        #print "So creating one!"
        newBCEntry = emptyBoundaryConditionCreate(bcEntry="velocityBCs")
    if "vcList" not in globalDict[bcEntry]:
        globalDict[bcEntry]["vcList"]=[]  # quietly create the necessary list

    # structure here is dict->list->dict->list
    vcListList = globalDict[bcEntry]["vcList"] # get list from BC dictionary (which is attached to main dictionary)

    wallVCdict = dict()
    wallVCdict["type"]="WallVC"
    wallVCdict["wall"]=wall
    wallVCdict["variables"]=[]


    vx=dict()
    vx["name"]="vx"
    vx["type"]="double"
    vx["value"]=0.0
    vy=dict()
    vy["name"]="vy"
    vy["type"]="double"
    vy["value"]=0.0
    vz=dict()
    vz["name"]="vz"
    vz["type"]="double"
    vz["value"]=0.0

    wallVCdict["variables"].append(vx)
    wallVCdict["variables"].append(vy)
    if dim == 3:
        wallVCdict["variables"].append(vz)

    vcListList.append(wallVCdict)

    return globalDict[bcEntry]


def wallSetFuncCreate(bcEntry="", wall="top", func="", **params):
    """
    Set a wall to be determined by a function for a Cartesian Mesh
    Args:
        func  : A function name to set the boundary condition
                One of:
                   [ Velocity_Lid_RampWithCentralMax, Velocity_SinusoidalLid ]

                These require the 'StgFEM_StandardConditionFunctions' plugin to be loaded
                and that the velocity FeVariable be named exactly 'VelocityField'
    """

    if func=="":
        _uw.utils.sendWarning("No function provided. A function must be set")

    globalDict = _uw.dictionary.GetDictionary()

    if bcEntry=="":
        field = globalDict["info"]["velocityField"]
        bcEntry = globalDict["info"][field+"BCs"]

    # Is it a good idea to make sure we have things like this on the dictionary always?
    # it makes some things a lot easier though
    if "dim" in globalDict.keys():
        dim = globalDict["dim"]
        #_uw.utils.sendWarning("Dim set to"+" "+str(dim))
    else:
        dim = 2    # If it's not set I am going to assume it's 2D
        _uw.utils.sendWarning("The value 'dim' not set in dictionary. So assuming dim=2")

    if bcEntry not in globalDict:
        #print "No boundary condition entry found with that name"
        #print "So creating one!"
        newBCEntry = emptyBoundaryConditionCreate(bcEntry=bcEntry)
    if "vcList" not in globalDict[bcEntry]:
        globalDict[bcEntry]["vcList"]=[]  # quietly create the necessary list

    # structure here is dict->list->dict->list
    vcListList = globalDict[bcEntry]["vcList"] # get list from BC dictionary (which is attached to main dictionary)

    wallVCdict = dict()
    wallVCdict["type"]="WallVC"
    wallVCdict["wall"]=wall
    wallVCdict["variables"]=[]

    # no slip by default
    vx=dict()
    vx["name"]="vx"
    vx["type"]="double"
    vx["value"]=0.0
    vy=dict()
    vy["name"]="vy"
    vy["type"]="double"
    vy["value"]=0.0
    vz=dict()
    vz["name"]="vz"
    vz["type"]="double"
    vz["value"]=0.0

    if wall=="left" or wall=="right":
        vy["type"]="func"
        vy["value"]=func
        wallVCdict["variables"].append(vx)
        wallVCdict["variables"].append(vy)
    if wall=="top" or wall=="bottom":
        vx["type"]="func"
        vx["value"]=func
        wallVCdict["variables"].append(vx)
        wallVCdict["variables"].append(vy)

    if dim == 3:
        if wall=="front" or wall=="back":
            vz["type"]="func"
            vz["value"]=func
            wallVCdict["variables"].append(vy)
            wallVCdict["variables"].append(vx)
        # takes care of all cases
        wallVCdict["variables"].append(vz)

    vcListList.append(wallVCdict)


    for key in params.keys():
        _globalDict[key]=params[key]

    return globalDict[bcEntry]


def wallTemperatureCreate(bcEntry="", wall="bottom", value=0.0):
    """
    Set a temperature on a given 'wall' on a Cartesian Temperature mesh
    """

    globalDict = _uw.dictionary.GetDictionary()

    if bcEntry=="":
        field = globalDict["info"]["temperatureField"]
        bcEntry = globalDict["info"][field+"BCs"]

    if bcEntry not in globalDict:
        #print "No boundary condition entry found with that name"
        #print "So creating one!"
        # Maybe be quiet for now?
        newBCEntry = emptyBoundaryConditionCreate(bcEntry)
    if "vcList" not in globalDict[bcEntry]:
        globalDict[bcEntry]["vcList"]=[]  # quietly create the necessary list

    # structure here is dict->list->dict->list
    vcListList = globalDict[bcEntry]["vcList"] # get list from BC dictionary (which is attached to main dictionary)

    wallVCdict = dict()
    wallVCdict["type"]="WallVC"
    wallVCdict["wall"]=wall
    wallVCdict["variables"]=[]


    temp=dict()
    temp["name"]="temperature"   # is the mesh variable name
    temp["type"]="double"
    temp["value"]=str(value)

    wallVCdict["variables"].append(temp)

    vcListList.append(wallVCdict)

    return globalDict[bcEntry]

def temperatureICSinusoidalCreate(
        icEntry="",
        TopLayerCoord=1.0,
        TopLayerBC=0,
        BottomLayerCoord=0.0,
        BottomLayerBC=1,
        PerturbationAmplitude=0.1,
        HorizontalWaveNumber=1,
        VerticalWaveNumber=1
):
    """
    Set initial sinusoidal temperature on Mesh
    """

    globalDict = _uw.dictionary.GetDictionary()

    if icEntry=="":
        field = globalDict["info"]["temperatureField"]
        icEntry = globalDict["info"][field+"ICs"]

    if icEntry not in globalDict:
        #print "No boundary condition entry found with that name"
        #print "So creating one!"
        # Maybe be quiet for now?
        newICEntry = emptyBoundaryConditionCreate(icEntry)
    if "vcList" not in globalDict[icEntry]:
        globalDict[icEntry]["vcList"]=[]  # quietly create the necessary list

    globalDict["SinusoidalTempIC_TopLayerCoord"]      = str(TopLayerCoord)
    globalDict["SinusoidalTempIC_TopLayerBC"]           = str(TopLayerBC)
    globalDict["SinusoidalTempIC_BottomLayerCoord"] = str(BottomLayerCoord)
    globalDict["SinusoidalTempIC_BottomLayerBC"]      = str(BottomLayerBC)
    globalDict["SinusoidalTempIC_PerturbationAmplitude"]    = str(PerturbationAmplitude)
    globalDict["SinusoidalTempIC_HorizontalWaveNumber"] = str(HorizontalWaveNumber)
    globalDict["SinusoidalTempIC_VerticalWaveNumber"]     = str(VerticalWaveNumber)

    # structure here is dict->list->dict->list
    vcListList = globalDict[icEntry]["vcList"] # get list from BC dictionary (which is attached to main dictionary)

    wallVCdict = dict()
    wallVCdict["type"]="AllNodesVC"
    wallVCdict["variables"]=[]


    temp=dict()
    temp["name"]="temperatureMeshVariable"
    temp["name"]="temperature"   # is the mesh variable name
    temp["type"]="func"
    temp["value"]="LinearWithSinusoidalPerturbation"

    wallVCdict["variables"].append(temp)

    vcListList.append(wallVCdict)

    return globalDict[icEntry]
