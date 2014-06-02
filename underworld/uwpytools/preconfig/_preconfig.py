# preconfig package - set up stuff like context, toolboxes and plugins
##############################################################################
# Using this for testing.
# Some functions are mirrored from _utils.py
# but don't need Underworld to be running.
##############################################################################

#import uwpytools as _uw


'''
For testing without Underworld running.
'''

_globalDict = None

import confutils as utils
import collections as _collections

def initDict():
    
    global _globalDict
    
    if _globalDict == None:
        _globalDict = GetTemplateEmptyDict()

    return _globalDict

def setParameters(**params):
    
    global _globalDict

    for key in params.keys():
        if key=="outputPath":
            _globalDict["checkpointWritePath"]=params[key]+"/Checkpoints"
        _globalDict[key]=params[key]
    
    return

def importToolBox(toolBox):
    """
    Known ToolBoxes:
                  Underworld
                  Solvers
                  gLucifer
    """
    global _globalDict

    #todo: check for known ToolBoxes
    _globalDict["import"].append(toolBox)

    return

def addPlugin(plugin,  context="context", **other ):
    """
    Known Plugins:
               StgFEM_FrequentOutput
               StgFEM_CPUTime
               StgFEM_StandardConditionFunctions
               Underworld_Vrms      : requires a GaussSwarm and VelocityField to be set
               Underworld_PressureCalibration
    """
    global _globalDict

    #todo: check for known plugins?

    pdict={}
    pdict["Type"]=plugin
    pdict["Context"]=context

    for key in other.keys():
        pdict[key]=other[key]

    _globalDict["plugins"].append(pdict)
        
    return

def initDefaultParameters():

    global _globalDict
    
    output="output"

    _globalDict["dim"]=2
    _globalDict["maxX"]=1.0
    _globalDict["maxY"]=1.0
    _globalDict["maxZ"]=1.0
    _globalDict["minX"]=0.0
    _globalDict["minY"]=0.0
    _globalDict["minZ"]=0.0
    _globalDict["resX"]=4
    _globalDict["resY"]=4
    _globalDict["resZ"]=1
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
    NewComponentEntryInStgDict( _globalDict,
                                name = "context",
                                Type = "UnderworldContext"
                                )
    _globalDict["penaltyNumber"]=0.1
    _globalDict["mgLevels"]=3
    _globalDict["saveDataEvery"]=1
    _globalDict["checkpointEvery"]=1
    _globalDict["checkpointWritePath"]=output+"/Checkpoints"
    _globalDict["checkpointAppendStep"]=1
    _globalDict["outputSlimmedXML"]=True

    # I think having these might make life easier for everyone - these are just names/flags
    _globalDict["context"]="context"
    _globalDict["gaussIntSwarm"]=0
    _globalDict["picIntSwarm"]=0
    _globalDict["FeVariableQ1"]=0 # maybe too abstract?
    _globalDict["FeVariableP0"]=0
    _globalDict["velocityFeVariable"]=0 # use these for now <-- will be set by mesh Creation 
    _globalDict["pressureFeVariable"]=0
    _globalDict["solver"]=0
    return

def GetTemplateEmptyDict():
    """
    Returns an empty python template dictionary for the user to fill as required.
    
    Args:
        None
    Returns:
        templateDict (OrderedDict): Python template ordered dictionary to utilise in StGermain.
    """
    dictionary = _collections.OrderedDict()
    dictionary["import"]=[]
    dictionary["plugins"]=[]
    dictionary["components"]=_collections.OrderedDict()
    return dictionary

def NewComponentEntryInStgDict( globalDict, **cptDefinitionArgs ):
    """
    Creates a new item in the components section of the dictionary

    Args:
        globalDict (dictionary) : The runtime global dictionary
        name=(string)            : The name of the component (must be unique) 
        Type=(string)             : The component type
        XXX=(some data)            : As defined by the component itself

    Returns:
        stgDictionary dictionary node for the new component 
    """


    if not "name" in cptDefinitionArgs:
        print "Error: NewComponentEntryInStgDict requires a (unique) 'name' "
        return
    
    if not "Type" in cptDefinitionArgs:
        print "Error: NewComponentEntryInStgDict requires a (valid) 'Type' "
        return


    globalDict["components"][cptDefinitionArgs["name"]] = cptDefinitionArgs


    return globalDict["components"][cptDefinitionArgs["name"]]


def GetCurrentPythonDictionary():
    """
    Returns the global copy of the Python equivalent of the StGermain dictionary. 

    Args:
        None
    Returns:
        dict (OrderedDict): Python ordered dictionary representation of the global StGermain dictionary.

    """

    global _globalDict
        return(_globalDict)

def PrettyDictionaryPrint(dict, indent=3):
    """
    Prints the provided dictionary in a human format.

    Args:
        dict (dict): The dictionary to print
        indent (unsigned): (optional) Required indent level
    Returns:
        PrettyDict (str):  The formated dictionary as a string.
    """

    import json
    print json.dumps(dict, indent=indent)
