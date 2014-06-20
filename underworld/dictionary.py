# helper tools for building an Underworld dictionary

import collections as _collections

def GetDictionary():
    """
    Returns the Underworld dictionary.

    Args:
        None
    Returns:
        dict (OrderedDict): Underworld dictionary.  This will be 'None' until underworld.Init() is called.

    """

    global _globalDict

    return _globalDict


def SetDictionary( pyDict ):
    """
    Sets / Overwrites the global copy of the python global dictionary.

    Args:
        dict (OrderedDict): Python ordered dictionary in the required format (be careful !)

    Returns:
        dict (OrderedDict): Returns the set dictionary
    """

    global _globalDict

    _globalDict = pyDict


def ClearDictionary():
    """
    Clears the global Underworld dictionary.

    Args:
        None

    Returns:
        Nothing
    """

    global _globalDict

    if _globalDict != None:
        _globalDict.clear()

    return

def GetConformingDictionary(uwdict={}):
    """
    Checks the provided dictionary has the required structure elements.
    Returns a new dictionary with all required structure (and content of provided dictionary).

    Args:
        dict (dict,OrderedDict) : Dictionary to check. Default is None.
    Returns:
        dict (OrderedDict)      : New conforming dictionary
    """
    if not isinstance(uwdict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    dictionary = _collections.OrderedDict(uwdict)
    if not "plugins" in dictionary:
        dictionary["plugins"]=[]
    if not "import" in dictionary:
        dictionary["import"]=[]
    if not "parameters" in dictionary:
        dictionary["parameters"]={}
    if not "info" in dictionary:
        dictionary["info"]={}

    if not "components" in dictionary:
        dictionary["components"]=_collections.OrderedDict()
    else:
        dictionary["components"]=_collections.OrderedDict(dictionary["components"])

    return dictionary


def UpdateDictWithComponent( uwdict=None, **cptDefinitionArgs ):

    """
    Updates the global python dictionary components with the information supplied.
    Will add a component if not found within dictionary.

    Args:
        uwdict (dict)            : Dictionary to update.  Default is globalDict.
        name=(string)            : The name of the component (must be unique)
        Type=(string)            : The component type
        mergeType (string)       : The equivalent of the U/W mergeType (replace or merge)
        XXX=(some data)          : As defined by the component itself

    Returns:
        stgDictionary dictionary node for the (new) component
    """

    if not uwdict:
        uwdict = GetDictionary()
    if not "name" in cptDefinitionArgs:
        print "Error: UpdateDictWithComponent requires a (unique) 'name' "
        print "Cpt: ", cptDefinitionArgs
        return

    if not "Type" in cptDefinitionArgs:
        print "Error: UpdateDictWithComponent requires a (valid) 'Type' "
        print "Cpt: ", cptDefinitionArgs
        return

    if "mergeType" in cptDefinitionArgs and cptDefinitionArgs["mergeType"] == "replace":
        del uwdict["components"][cptDefinitionArgs["name"]]
        uwdict["components"][cptDefinitionArgs["name"]] = cptDefinitionArgs
    else:
        if cptDefinitionArgs["name"] in uwdict["components"]:
            uwdict["components"][cptDefinitionArgs["name"]].update(cptDefinitionArgs)
        else:
            uwdict["components"][cptDefinitionArgs["name"]] = cptDefinitionArgs


    return uwdict["components"][cptDefinitionArgs["name"]]



def setParameters(**params):

    global _globalDict

    ## No way this special case can be left in here !!!!!

    for key in params.keys():
        if key=="outputPath":
            _globalDict["checkpointWritePath"]=params[key]+"/Checkpoints"
            _globalDict["checkpointReadPath"]=params[key]+"/Checkpoints"
        _globalDict[key]=params[key]

    return

def importToolBox(toolBox):
    """
    Known ToolBoxes:
                  Underworld
                  Solvers     (this is supposed to be the default in 2.0)
                  gLucifer
                  viscoelastic (this is supposed to be the default in 2.0)
    """

    global _globalDict

    #todo: check for known ToolBoxes
    #todo: check config to see which toolboxes are available
    _globalDict["import"].append(toolBox)

    return

def getInfo():

    return GetDictionary()["info"]

def addPlugin(plugin,  context="context", **other ):
    """
    Known Plugins:
               StgFEM_FrequentOutput                  - (deprecated)
               StgFEM_CPUTime                         - (deprecated)
               StgFEM_StandardConditionFunctions
               Underworld_Vrms      : requires a GaussSwarm and VelocityField to be set - (deprecated)
               Underworld_PressureCalibration
    """
    global _globalDict

    #todo: check for known plugins?
    #todo: eliminate plugins completely.

    pdict={}
    pdict["Type"]=plugin
    pdict["Context"]=context

    for key in other.keys():
        pdict[key]=other[key]

    _globalDict["plugins"].append(pdict)

    return

# These should not be set here, they should be in the initialisation for each of the modules
# and they should be in the parameters list of the dictionary
# e.g.
# globalDict["parameters"]["geometry"]["dim"] = 2
def addCheckPointVariables(checkVars=[]):
    """
    Add FieldVariables to Checkpoint
    """
    global _globalDict

    # For the uw parameters like "outputPath", they should also be moved to a global list
    if "FieldVariablesToCheckpoint" in _globalDict:
        for item in checkVars:
            _globalDict["FieldVariablesToCheckpoint"].append(item)
    else:
        _globalDict["FieldVariablesToCheckpoint"]=checkVars

    return

# We want o put these somewhere else..
# But for the moment I think UW requires them at top level in dict.
def addSaveVariables(saveVars=[]):
    """
    Add FieldVariables to Save
    """
    global _globalDict

    if "FieldVariablesToSave" in _globalDict:
        for item in saveVars:
            _globalDict["FieldVariablesToSave"].append(item)
    else:
        _globalDict["FieldVariablesToSave"]=saveVars

    return


def initDefaultParameters():

    global _globalDict

    output="output"

 #   _globalDict["dim"]=2
 #  _globalDict["maxX"]=1.0
  #  _globalDict["maxY"]=1.0
  #  _globalDict["maxZ"]=1.0
   # _globalDict["minX"]=0.0
  #  _globalDict["minY"]=0.0
   # _globalDict["minZ"]=0.0
   # _globalDict["resX"]=4
   # _globalDict["resY"]=4
   # _globalDict["resZ"]=1
    _globalDict["outputPath"]=output
    _globalDict["shadowDepth"]=1
    _globalDict["maxTimeSteps"]=1
    _globalDict["courantFactor"]=0.25
    _globalDict["particlesPerCell"]=20
    _globalDict["seed"]=13
    _globalDict["allowUnbalancing"]=True
    _globalDict["buildElementNodeTbl"]=True
    _globalDict["particleLayoutType"]="random"
    _globalDict["dumpEvery"]= 1
    _globalDict["import"].append("Underworld")
    UpdateDictWithComponent(    name = "context",
                                Type = "UnderworldContext"
                                )
    _globalDict["penaltyNumber"]=0.1
    _globalDict["mgLevels"]=3
    _globalDict["saveDataEvery"]=1
    _globalDict["checkpointEvery"]=1
    _globalDict["checkpointWritePath"]=output+"/Checkpoints"
    _globalDict["checkpointReadPath"]=output+"/Checkpoints"
    _globalDict["checkpointAppendStep"]=0
    _globalDict["outputSlimmedXML"]=True

    _globalDict["gaussParticlesX"]=2
    _globalDict["gaussParticlesY"]=2
    _globalDict["gaussParticlesZ"]=2

    # I think having these might make life easier for everyone - these are just names/flags
    _globalDict["context"]="context"
    _globalDict["gaussIntSwarm"]=0
    _globalDict["picIntSwarm"]=0
    _globalDict["FeVariableQ1"]=0 # maybe too abstract?
    _globalDict["FeVariableP0"]=0
    _globalDict["velocityFeVariable"]=0 # use these for now <-- will be set by mesh Creation
    _globalDict["pressureFeVariable"]=0
    _globalDict["solver"]=0

    addPlugin("StgFEM_FrequentOutput")
    addPlugin("StgFEM_CPUTime")
    addPlugin("StgFEM_StandardConditionFunctions")

    return

def PrintPretty(dict=None, indent=3):
    """
       Prints the provided dictionary in a human format.

       Args:
       dict (dict): The dictionary to print. Default is the global dictionary.
       indent (unsigned): (optional) Required indent level
       Returns:
       PrettyDict (str):  The formated dictionary as a string.
       """

    import json
    if not dict:
        global _globalDict
        dict = _globalDict
    print json.dumps(dict, indent=indent)


_globalDict = None
