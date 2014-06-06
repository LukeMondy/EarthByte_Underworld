# Solvers package - activates and configures the "new solver system"

import underworld
import os
import glob

'''
Um ...

'''

def pathToBuild():

    # If the build is separated from the source, then usually this will need to be set
    uwpathvariable = os.environ["UW_DIR"]  
    uwexecutable = os.path.join(uwpathvariable,"bin","Underworld")

    # I think this is about how much effort we should put in - somebody set this 
    # variable and it point to an executable file ...  

    if os.path.exists(uwpathvariable) and os.path.exists(uwexecutable) and os.access(uwexecutable, os.X_OK):
        return uwpathvariable

    # Otherwise, construct a path relative to this module file     (assuming in underworld/utils)

    # Should use the path-join routines - this is Unix specific

    uwrootpathvariable = os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))

    # get / count the number of valid installations
    uwconfigpathlist = glob.glob(os.path.join(uwrootpathvariable,"*build*","config.cfg"))

    if len(uwconfigpathlist) <> 1:
        print "Unable to find a unique build path - try setting the $UW_DIR environment variable"
        return ""

    uwpathvariable = os.path.dirname(uwconfigpathlist[0])    
    uwexecutable = uwpathvariable + "/bin/Underworld"

    if os.path.exists(uwexecutable) and os.access(uwexecutable, os.X_OK):
        return uwpathvariable

    return ""    


def pathToConfigFile():

    uwpathvariable = pathToBuild()

    # Check for invalid result - if triggered, then return too
    if (uwpathvariable == ""):
        return ""

    configfile = os.path.join(uwpathvariable,"config.cfg")

    if (os.path.exists(configfile)):
        return configfile
    else:
        print "Unable to locate config.cfg file - does $UW_DIR point to a valid build ?"
        return ""


def configDictFromConfigFile():
    
    configFileName = pathToConfigFile()
    
    if (configFileName == ""):
        print "Unable to locate config.cfg file"
        return ""

    configFile = open(configFileName,"r")
    rawconfig = configFile.read()    
    lineslist = rawconfig.splitlines()

    confDict={}
    for string in lineslist:
        confKey,separator,confValue = string.partition(' = ')
        confDict[confKey] = confValue
    
    return confDict

# This is not just a check !!
def checkForNewComponentName(globalDict, componentName):

    """
    Checks if a component is already in the dictionary.
    If it is then returns a new name else returns same name.
    """

    if componentName in globalDict["components"]: # We could just make a new name or refuse to add here.
        # I might try making a new name
        #print sendWarning("A component named "+componentName+"  already exists. ")
        count=1
        while componentName in globalDict["components"]:
            count += 1
            componentName = componentName+str(count)
        #print sendWarning("Returning name:  "+componentName)
    return componentName



## Should make a point of using globalDict for this 

def uniqueComponentNameGlobalDict(componentName):

    globalDict = underworld.GetCurrentPythonDictionary()

    if componentName in globalDict["components"]:             # We could just make a new name or refuse to add here.
        count=1

        while componentName in globalDict["components"]:
            count += 1
            componentName = componentName+str(count)

    return componentName


def warnMissingComponent(globalDict, componentName):
    """
    Checks if a component is in the dictionary.
    If it is not, then prints a warning
    """
    if componentName not in globalDict["components"]:
        print sendWarning("The "+componentName+" is not in the dictionary.")
        return 1
    return 0

def sendWarning(message):
    green="\033[0;32m"
    endcol="\033[00m"
    boldgreen="\033[1;32m"
    print " "+green+"*  Warning: "+endcol+boldgreen+message+endcol
    return

def sendError(message):
    purple="\033[0;35m"
    endcol="\033[00m"
    boldpurple="\033[1;35m"
    print " "+boldpurple+"*  Error: "+message+endcol
    return

def listComponentsByType(globalDict, Type):
    
    compList=[]
    comps=globalDict["components"]
    for compkey in comps:
        if Type in comps[compkey]["Type"]:  # looser test than ==
            compList.append(comps[compkey])
    
    return compList
