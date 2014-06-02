# preconfig package - set up stuff like context, toolboxes and plugins
# Using this for testing.
# Some functions are mirrored from _utils.py
# but don't need Underworld to be running.
import collections as _collections
#import utils as _utils

#_globalDict = None
##############################################################################
## This code adds what is required to the python dictionary 
## to set up a Mesh for Underworld.
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary for a Mesh
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''
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

def warnMissingComponent(globalDict, componentName):
    """
    Checks if a component is in the dictionary.
    If it is not, then prints a warning
    """
    if componentName not in globalDict["components"]:
        #print "\033[0;32mWarning\033[00m\033[1;32m: The "+componentName+" is not in the dictionary.\033[00m"
        print sendWarning("The "+componentName+" is not in the dictionary.")
        return 1
    return 0

def sendWarning(message):
    green="\033[0;32m"
    endcol="\033[00m"
    boldgreen="\033[1;32m"
    print green+"Warning: "+endcol+boldgreen+message+endcol
    return

def sendError(message):
    purple="\033[0;35m"
    endcol="\033[00m"
    boldpurple="\033[1;35m"
    print boldpurple+"Error: "+message+endcol
    return

def listComponentsByType(globalDict, Type):
    
    compList=[]
    comps=globalDict["components"]
    for compkey in comps:
        if Type in comps[compkey]["Type"]:  # looser test than ==
            compList.append(comps[compkey])
    
    return compList

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


def checkForNewComponentName(globalDict, componentName):
    """
    Checks if a component is already in the dictionary.
    If it is then returns a new name else returns same name.
    """
    if componentName in globalDict["components"]: # We could just make a new name or refuse to add here.
        # I might try making a new name
        #print "\033[0;32mWarning\033[00m\033[1;32m: A component of that name already exists. Making a new one with a different name\033[00m"
        print sendWarning("A component of that name already exists. Making a new one with a different name")
        count=1
        while componentName in globalDict["components"]:
            count += 1
            componentName = componentName+str(count)

    return componentName
