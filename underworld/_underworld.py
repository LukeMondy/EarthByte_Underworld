import sys as _sys
import os as _os
import tempfile as _tempfile
import subprocess as _subprocess
import xml.etree.cElementTree as _ET
import collections as _collections
from libUnderworld import *
import utils
import fields

def _elemGetKey(elem):
    if 'name' in elem.attrib:
        return elem.attrib['name']
    else:
        return elem.tag.split("}")[1]

def _elementToDict(elem):
    # set key
    dicttag = elem.tag.split("}")[1]
    # get children count
    if(len(elem)):
        value = _collections.OrderedDict()
        for child in elem:
            if child.tag.split('}')[1] in ['param', 'list', 'struct']:  # add this check to ensure we use only the valid guys
                childkey = _elemGetKey(child)
                valueGuy = _elementToDict(child)
                if dicttag in ['struct','StGermainData']:
                    if not childkey in value:
                        # if no entry with this key, add one
                        value[childkey] = valueGuy
                elif dicttag == 'list' :
                    # if entry already exists,
                    try:
                        # assume it's a list, and add item to list
                        value[childkey].append(valueGuy)
                    except:
                        # if that fails, lets create a list and add item
                        value[childkey] = []
                        value[childkey].append(valueGuy)
                else: 
                    print "Error.. param cannot contain sub values"
                    print childkey, valueguy
                    _sys.exit(2)
        # now, if elem is a list, ensure there's only one dict entry, and then simply return the list
        if dicttag == "list":
            if len(value) == 1:
                tempval = value.items()[0][1]
                value = tempval
            else:
                print "Error: list {} contains more than one item type. Ignoring.".format(_elemGetKey(child))

    else:
        value = _convToNativeWherePossible(elem.text)

    return value

def _convToNativeWherePossible(string):
    if string == None:
        return None
    elif string.lower()=="true":
        return True
    elif string.lower()=="false":
        return False
    else:
        try:
            return int(string)
        except ValueError:
            try:
                return float(string)
            except ValueError:
                return string


def _dictToUWElementTree(inputDict):
    # lets create root element
    root = _ET.Element('StGermainData')
    root.attrib['xmlns'] = 'http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003'

    # now to add subElements
    for k,v in inputDict.items():
        _itemToElement(v,k,root)

    return root


def _itemToElement(inputItem,inputItemName,inputEl):
    if type(inputItem) == list:
        subEl = _ET.SubElement(inputEl, 'list')
        if inputItemName != '':
            subEl.attrib['name']=inputItemName
        for item in inputItem:
            _itemToElement(item,'',subEl)
    elif type(inputItem) in [_collections.OrderedDict,dict]:
        subEl = _ET.SubElement(inputEl, 'struct')
        if inputItemName != '':
            subEl.attrib['name']=inputItemName
        for k,v in inputItem.items():
            _itemToElement(v,k,subEl)
    elif type(inputItem) in [str,float,int,bool,unicode]:
        subEl = _ET.SubElement(inputEl, 'param')
        if inputItemName != '':
            subEl.attrib['name']=inputItemName
        subEl.text = str(inputItem)
    elif not inputItem:
        subEl = _ET.SubElement(inputEl, 'param')
        if inputItemName != '':
            subEl.attrib['name']=inputItemName
        subEl.text = "\t"        
    else:
        print "Error.. Unknown type encountered"
        print "key =", inputItemName
        print "value =", inputItem
        _sys.exit(2)


def _indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def GetCurrentDictionary():
    """
    Gets a copy of the current StGermain dictionary
    This routines can only be utilised after the Init phase has been completed.

    Args:
        None
    Returns:
        dict (OrderedDict): Python ordered dictionary representation of StGermain dictionary.

    """
    global _data
    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    xmlString = StGermain_Tools.StgDictAsXMLString( _data )
    root = _ET.fromstring(xmlString)
    dict = _elementToDict(root)
    if not dict:
        dict = GetTemplateEmptyDict()

    # lets remove sources, as this information is lost anyhow.  if there's a need, we'll come back to it. 
    if "sources" in dict:
        del dict["sources"]
    return dict


### These functions allow us to build, merge and manipulate the global dictionary as a python object
### and, later, it can be converted back to XML strings for the blessed StGermain to injest.     


def GetCurrentPythonDictionary():
    """
    Returns the global copy of the Python equivalent of the StGermain dictionary. If the dictionary
    does not have a "parameters" or "info" node then this is added

    Args:
        None
    Returns:
        dict (OrderedDict): Python ordered dictionary representation of the global StGermain dictionary.

    """

    global _globalDict

    # 
    # This should never be the case since _globalDict is set as soon as an XML dictionary is created by the Init()
    # 

    if _globalDict == None:
        _globalDict = GetCurrentDictionary()
        if "parameters" not in _globalDict.keys():
            _globalDict["parameters"] = {}
        if "info" not in _globalDict.keys():
            _globalDict["info"] = {}

    return(_globalDict)            


def SetCurrentPythonDictionary( pyDict ):
    """
    Sets / Overwrites the global copy of the python global dictionary. 

    Args:
        dict (OrderedDict): Python ordered dictionary in the required format (be careful !)

    Returns:
        None
    """

    global _globalDict

    _globalDict = pyDict

    return            


def PushCurrentPythonDictionaryToStg():
    """
        Loads the current global Python dictionary into the XML dictionary of StGermain

    """

    global _globalDict

    SetDictionary( _globalDict )

    return



## We might want to move these basic functions to underworld.utils so we are not importing so much

## Also, I'm not a fan of these multi-type functions. 
## I think it makes everything a lot less clear !

def StgXMLFileAsPyDictionary( xmlFile, flatten=False ):
    """
    Converts the provided XML file into a python dictionary.

    Args:
        xmlFile (str, File): File to convert.  Either an opened file, or the filename as a string. 
        flatten (Bool): if xmlFile is a string, the flattenXML command may be called to flatten the file. 
    Returns:
        dict (OrderedDict): Python ordered dictionary containing the file contents.
    
    """

    if flatten:
        if (type(xmlFile) == file):
            print "Unable to flatten a file object - this option requires the file path"        
        else:
            # Use the StGermain FlattenXML executable 
            build_path = utils.pathToBuild()
            flattenXML = _os.path.join(build_path,"bin","FlattenXML")
            current_dir = _os.getcwd()
            tempdir = _tempfile.mkdtemp()
            _os.chdir(tempdir)
            _subprocess.call([flattenXML,xmlFile])
            xmlFile = _os.path.join(tempdir,"output.xml")
            _os.chdir(current_dir)

    if type(xmlFile) == str:
        theFile = open(xmlFile,'r') 
    elif type(xmlFile) == file:
        theFile = xmlFile
    else:
        print "You must pass in a file or a filename"
        return

    xmlString = theFile.read().replace('\n', '')    
    root = _ET.fromstring(xmlString)
    dict = _elementToDict(root)
    if type(xmlFile) == str:
        theFile.close()


    ## remove any temp files ... 
    ## HERE

    return dict

def WritePyDictToJSONFile( theDict, jsonFile ):
    """
    Converts the provided python dictionary into a JSON file.
    If running in parallel, make sure only one processor executes this command. 

    Args:
        theDict (dict): Python dictionary to convert.
        jsonFile (str, File): File to write to.  Either an opened file, or the filename as a string. 
    
    """
    import json
    if type(jsonFile) == str:
        theFile = open(jsonFile,'w') 
    elif type(jsonFile) == file:
        theFile = jsonFile
    else:
        print "You must pass in a file or a filename"
        return
    theFile.write(json.dumps(theDict, indent=3))
    # close the file if we opened it
    if type(jsonFile) == str:
        theFile.close()

def ReadJSONFileToPyDict( jsonFile ):
    """
    Converts the provided JSON file into a python dictionary

    Args:
        jsonFile (str, File): File to convert.  Either an opened file, or the filename as a string. 
    
    Returns:
        dict (OrderedDict): Python ordered dictionary containing the file contents.
    """

    import json
    if type(jsonFile) == str:
        theFile = open(jsonFile,'r') 
    elif type(jsonFile) == file:
        theFile = jsonFile
    else:
        print "You must pass in a file or a filename"
        return
    theDict = json.load(theFile, object_pairs_hook=_collections.OrderedDict)
    # close the file if we opened it
    if type(jsonFile) == str:
        theFile.close()
    return theDict

def GetTemplateEmptyDict():
    """
    Returns an empty python template dictionary for the user to fill as required.
    
    Args:
        None
    Returns:
        templateDict (OrderedDict): Python template ordered dictionary to utilise in StGermain.
    """
    dictionary = _collections.OrderedDict()
    dictionary["plugins"]=[]
    dictionary["import"]=[]
    dictionary["components"]=_collections.OrderedDict()
    return dictionary

def SetDictionary( pyDict, tag=_os.path.realpath(_sys.argv[0]) ):
    """
    Sets the provided python dictionary as the StGermain dictionary.
    This routines can only be utilised after the Init phase has been completed, and 
    should usually be run before the Construct/Build/Initialise phases.

    Args:
        pyDict (dict): Python dictionary to set as StGermain dictionary.
        tag (str): (optional) Tag to add to StGermain sources dictionary for this object.  Defaults to the script name.
    """
    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == True or _isBuilt == True :
        print "StGermain has already been Constructed/Built/Initialised."
        print "Modifications to the dictionary will possibly have no effect."
        print "You may need to restart, either by quitting or calling Finalise()"

    root = _dictToUWElementTree(pyDict)
    xmlString = _ET.tostring(root,encoding = 'utf-8', method = 'xml')
    StGermain_Tools.StgSetDictFromXMLString(_data,xmlString, tag)



##
##    

def AddToDictionary( pyDict, tag=_os.path.realpath(_sys.argv[0]) ):
    """
    Adds the provided python dictionary to the StGermain dictionary.
    This routines can only be utilised after the Init phase has been completed, and 
    should usually be run before the Construct/Build/Initialise phases.

    Args:
        pyDict (dict): Python dictionary to add to StGermain dictionary.
        tag (str): (optional) Tag to add to StGermain sources dictionary for this object.  Defaults to the script name.
    """

    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == True or _isBuilt == True:
        print "StGermain has already been Constructed/Built/Initialised."
        print "Modifications to the dictionary will possibly have no effect."
        print "You may need to restart, either by quitting or calling Finalise()"

    root = _dictToUWElementTree(pyDict)
    xmlString = _ET.tostring(root, encoding = 'utf-8', method = 'xml')
    StGermain_Tools.StgAddToDictFromXMLString(_data, xmlString, tag)

    return

def InitFromCommandLine(extraArgs=[]):
    """
    Runs the StGermain init phase using the arguments provided at the script command line.

    Args:
        extraArgs (list, str): (optional) List of extra arguments, provided as a list or as a string
    Returns:
        Nothing
    """

    global _data
    global _globalDict


    if _data != None:
        print "StGermain has already been initialised."
        print "Will finalise previous instance, and then re-init"
        Finalise()
        # And eliminate any stored dictionary to match the state of the underworld XML dictionary
        globalDict = GetCurrentPythonDictionary()
        globalDict.clear()

    if type(extraArgs) == str:
        extraArgs = extraArgs.split()

    InitWithArgs( _sys.argv[1:]+extraArgs )

    _globalDict = GetCurrentDictionary()
    if "parameters" not in _globalDict.keys():
        _globalDict["parameters"] = {}
    if "info" not in _globalDict.keys():
        _globalDict["info"] = {}

    return



def Init():
    """
    Runs the StGermain init phase without any arguments, initialises the global dictionary

    Args:
        None
    Returns:
        Nothing
    """

    global _data
    global _globalDict
 
    if _data != None:
        print "StGermain has already been initialised."
        print "Will finalise previous instance, and then re-init"
        Finalise()

    _data  = StGermain_Tools.StgInit( [_sys.argv[0]] )

    _globalDict = GetCurrentDictionary()
    if "parameters" not in _globalDict.keys():
        _globalDict["parameters"] = {}
    if "info" not in _globalDict.keys():
        _globalDict["info"] = {}

    return


def InitWithArgs( args ):
    """
    Runs the StGermain init phase using the arguments provided.  Will also load any JSON files.
    Note that any JSON files will be loaded last, and overwrite any previous items. 

    Args:
        args (list, str): List of arguments, provided as a list or as a string

    Returns:
        Nothing

    """

    global _data
    global _globalDict

    if _data != None:
        print "StGermain has already been initialised."
        print "Will finalise previous instance, and then re-init"
        Finalise()


    # if type is string, convert to list
    if type(args) == str:
        args = args.split()
    args.insert( 0, _sys.argv[0])  #insert executable name infront of user supplied args
    _data  = StGermain_Tools.StgInit( args )

    # lets add any json files now 

    for arg in args:
        if arg.endswith(('.json','.JSON')):
            jdict = ReadJSONFileToPyDict(arg)
            AddToDictionary(jdict, tag=arg)

    _globalDict = GetCurrentDictionary()
    if "parameters" not in _globalDict.keys():
        _globalDict["parameters"] = {}
    if "info" not in _globalDict.keys():
        _globalDict["info"] = {}



    return

def Construct( PushPythonDictionary=True ):
    """
    Runs the StGermain Construct simulation phase.
    An Init function must be called before this.

    The global Python dictionary is pushed by default since 
    this is the mechanism for carrying changes to the configuration.
    It would be better to enforce this by keeping the XML change functions
    out of the way or only using those. 

    Args:
        None
    Returns:
        Nothing
    """
    global _data
    global _isConstructed  # private

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == True:
        print "StGermain has already been constructed."
        print "To re-construct, you will need to restart either by quitting or calling Finalise()"
        return

    # Check the global python dictionary has been set - if it has not 
    # been retrieved at all then the first call to fetch it will populate it 
    # with the XML dictionary - so this test is self-determined !    

    globalDict = GetCurrentPythonDictionary()

    if PushPythonDictionary and globalDict != None:
        PushCurrentPythonDictionaryToStg()

    StGermain_Tools.StgConstruct(_data)
    _isConstructed = True

def BuildAndInitialise():
    """
    Runs the StGermain Build and Initialise phases for the simulation.
    An Init function must be called before this.

    Args:
        None
    Returns:
        Nothing
    """

    global _data
    global _isConstructed
    global _isBuilt # private
    global isBuilt  # public

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == False:
        print "StGermain has not been constructed. You will need to construct using Construct()"
        return
    if _isBuilt == True:
        print "StGermain has already been built."
        print "To re-build, you will need to restart either by quitting or calling Finalise()"
        return
    StGermain_Tools.StgBuildAndInitialise(_data)
    _isBuilt = True


def RunMainLoop():
    """
    Runs the StGermain main simulation loop.
    An Init function must be called before this, as well as the Construct, Build & Initialise routines.

    Args:
        None
    Returns:
        Nothing
    """

    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isBuilt == False:
        print "StGermain has not been built. You will need to build using BuildAndInitialise()"
        return
    StGermain_Tools.StgMainLoop(_data)

def Step(context=None, steps=1):
    """
    Performs required number of steps of the simulation. 
    An Init function must be called before this, as well as the Construct, Build & Initialise routines.
    This step will ignore all simulation flow control parameters (such as maxTimeStep).

    Args:
        context (Swig Context*): (optional) Context for current simulation.  
            If none provided, will search LiveComponentRegister    for component with name "context".
        steps (int): is the number of steps to take
        [TBA] time (float): is the maximum time to step forward 
            If both are supplied, the first termination condition met is the one acted upon.
    Returns:
        Nothing
    """
    global _data
    global _isConstructed
    global _isBuilt

    if context == None:
        context = GetLiveComponent("context")

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isBuilt == False:
        print "StGermain has not been built. You will need to build using BuildAndInitialise()"
        return

    for i in range(0,steps):
        if( context.needUpdate == 1 ):
            StGermain.AbstractContext_Update(context)

        StGermain.AbstractContext_Step(context, context.dt)


def RunAll():
    """
    Runs the entire StGermain simulation (Construct, Build, Initialise, Execute).
    An Init function must be called before this.

    Args:
        None
    Returns:
        Nothing
    """

    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isBuilt == True:
        print "StGermain has already been built."
        print "To re-build, you will need to restart either by quitting or calling Finalise()"
        return
    StGermain_Tools.StgRun(_data)

def Finalise():
    """
    Finalises / tears down the StGermain simulation.

    Args:
        None
    Returns:
        Nothing
    """
    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    StGermain_Tools.StgFinalise( _data )
    _resetAll()

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

'''  # currently not implemented
def SetOutputPath( outputPath ):
    global _data
    global _isConstructed
    global _isBuilt

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    StGermain_Tools.StgSetOutputPath( _data, outputPath )
'''

def rank():
    """
    Returns the rank of the current processors.

    Args:
        None
    Returns:
        rank (unsigned) : Rank of current processor.
    """
    global _data

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    return _data.rank

def nProcs():
    """
    Returns the number of processors being utilised by the simulation.

    Args:
        None
    Returns:
        nProcs (unsigned) : Number of processors.
    """

    global _data

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    return _data.nProcs

def _resetAll():
    global _data
    global _isConstructed
    global _isBuilt
    global _globalDict

    _data = None
    _isConstructed = False
    _isBuilt = False
    _globalDict.clear()

def GetLiveComponent( componentName ):
    """
    Returns a live component from the StGermain live component register.
    Note: This can only be run after the Construct phase.
    
    Args:
        componentName (str):  Name of the component you are trying to retrieve
    Returns:
        component (swig*) : A swig pointer type to the requested component, or None if not found.
    """

    global _data
    global _isContructed
    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == False:
        print "StGermain has not been constructed. You will need to construct using Construct()"
        return
    try:
        return StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), componentName )
    except:
        print "Component \'%s\' not found in the live component register." % componentName
        return None

def GetLiveComponentsAsDictionary():
    """
    Returns a dictionary of ALL live components organised by name [as component (swig*) ]
    """

    global _data
    global _isContructed
    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return
    if _isConstructed == False:
        print "StGermain has not been constructed. You will need to construct using Construct()"
        return

    thisDict = GetCurrentDictionary()

    componentsDict = {}
    for cptName in thisDict["components"]:

        try: 
            liveCpt = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), cptName )
            componentsDict[cptName] = liveCpt

        except:
            liveCpt = 0 

    return componentsDict


def FieldVariable_InterpolateValueAt( fieldVar, coord ):
    """
    Returns the interpolated result of a field variable.
    
    Args:
        fieldVar (str, Swig FieldVariable*): Field variable to query, either provided by name (str) or with a swig pointer.
        coord (tuple(float)): The location to query
    Returns:
        result (list(float)): The interpolated result as a list
        success (str): "OTHER_PROC"     - Query location is on another processor
                       "LOCAL"          - Query location is on this processor
                       "SHADOW"         - Query location is in shadow space of this processor
                       "OUTSIDE_GLOBAL" - Query location is outside the global domain of this field variable
    """
    import c_arrays
    
    if type(fieldVar)==str:
        fieldVar = GetLiveComponent(fieldVar)

    result = c_arrays.DoubleArray(fieldVar.fieldComponentCount)

    InterpolationResult = StgDomain.FieldVariable_InterpolateValueAt( fieldVar, coord, result.cast() )

    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"

    toreturn = []
    for i in range(0,fieldVar.fieldComponentCount):
        toreturn.append(result.__getitem__(i))

    return toreturn, InterpolationResult


# This is the "proper" way to do things ... differentiate the kernel etc etc.
# But this is not particularly stable and needs a lot of points to work.
# We probably need a more sophisticated weighting scheme to make this work
# But calculating a FD gradient is quite accurate ... 

def RBFGradientVariable_InterpolateValueAt1( RBFfieldVar, coord , axis ):
    """
    Returns the interpolated result of an RBF variable.
    
    Args:
        fieldVar (str, Swig FieldVariable*): Field variable to query, either provided by name (str) or with a swig pointer.
        coord (tuple(float)): The location to query
    Returns:
        result (list(float)): The interpolated result as a list
        success (str): "OTHER_PROC"     - Query location is on another processor
                       "LOCAL"          - Query location is on this processor
                       "SHADOW"         - Query location is in shadow space of this processor
                       "OUTSIDE_GLOBAL" - Query location is outside the global domain of this field variable
    """
    
    import c_arrays
    
    if type(RBFfieldVar)==str:
        RBFfieldVar = GetLiveComponent(RBFfieldVar)

    result = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    coordArray = c_arrays.DoubleArray(3)

    coordArray[0] = coord[0]
    coordArray[1] = coord[1]
    coordArray[2] = coord[2]

    InterpolationResult = Underworld.RBFFieldVariable_InterpolateGradientValueAt( RBFfieldVar, coordArray.cast(), result.cast(), axis )

    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"

    toreturn = []
    for i in range(0,RBFfieldVar.fieldComponentCount):
        toreturn.append(result.__getitem__(i))

    return toreturn, InterpolationResult


# 
# This is a Finite Difference calculation of the gradient of the RBF field 

def RBFGradientVariable_InterpolateValueAt( RBFfieldVar, coord , axis ):
    """
    Returns the interpolated result of an RBF variable.
    
    Args:
        fieldVar (str, Swig FieldVariable*): Field variable to query, either provided by name (str) or with a swig pointer.
        coord (tuple(float)): The location to query
    Returns:
        result (list(float)): The interpolated result as a list
        success (str): "OTHER_PROC"     - Query location is on another processor
                       "LOCAL"          - Query location is on this processor
                       "SHADOW"         - Query location is in shadow space of this processor
                       "OUTSIDE_GLOBAL" - Query location is outside the global domain of this field variable
    """
    
    import c_arrays
    
    if type(RBFfieldVar)==str:
        RBFfieldVar = GetLiveComponent(RBFfieldVar)

    result0 = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    resultM = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)
    resultP = c_arrays.DoubleArray(RBFfieldVar.fieldComponentCount)

    coordArray = c_arrays.DoubleArray(3)
    coordArray[0] = coord[0]
    coordArray[1] = coord[1]
    # if()
    coordArray[2] = coord[2]

    # Compare it with this result

    supportSize = RBFfieldVar.rbfManager.particleSupportRadius
    
    InterpolationResult = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), result0.cast())

    if InterpolationResult == 3:
        return [], "OUTSIDE_GLOBAL"

    coordArray[axis] -= supportSize * 0.5
    InterpolationResultM = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), resultM.cast())

    coordArray[axis] += supportSize
    InterpolationResultP = Underworld._RBFFieldVariable_InterpolateValueAt( RBFfieldVar, coordArray.cast(), resultP.cast())

    
    # Weights 

    WM = -1.0; W0 =  0.0; WP =  1.0
    delta = supportSize

    if InterpolationResultM == 3:
        WM = 0.0
        W0 = -1.0
        delta = 0.5 * supportSize

    if InterpolationResultP == 3:
        WP = 0.0
        W0 += 1.0  # If both the sample points are outside the domain, the result is zero ... 
        delta = 0.5 * supportSize

    # if(axis == 1):
    #     print "Coord P - {} - {}".format(coordArray[1],coord[1])
    #     print "Weights {}/{}/{} Del {}".format(WM,W0,WP,delta)
    #     print "Values  {}/{}/{} ".format(resultM[0],result0[0],resultP[0])
    #     print "Results {}/{}/{} ".format(InterpolationResultM,InterpolationResult,InterpolationResultP)


    toreturn = []
    for i in range(0,RBFfieldVar.fieldComponentCount):
        toreturn.append( (WM * resultM[i] + W0 * result0[i] + WP * resultP[i]) / delta )


    if InterpolationResult == 0:
        InterpolationResult = "OTHER_PROC"
    elif InterpolationResult == 1:
        InterpolationResult = "LOCAL"
    elif InterpolationResult == 2:
        InterpolationResult = "SHADOW"
    elif InterpolationResult == 3:
        InterpolationResult = "OUTSIDE_GLOBAL"


    return toreturn, InterpolationResult



def FieldVariable_GetMinFieldMagnitude( fieldVar ):
    """
    Returns the minimum value of a field variable
    
    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (float): the minimum field variable magnitude
    """
    return StgDomain.FieldVariable_GetMinGlobalFieldMagnitude( fieldVar )

def FieldVariable_GetMaxFieldMagnitude( fieldVar ):
    """
    Returns the maximum value of a field variable
    
    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        max (float): the maximum field variable magnitude
    """
    return StgDomain.FieldVariable_GetMaxGlobalFieldMagnitude( fieldVar )

def FieldVariable_GetMinAndMaxLocalCoords( fieldVar ):
    """
    Returns the domain of this field variable local to this process as two tuples (min,max)
    
    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (tuple(float)): the minimum local coordinate range
        max (tuple(float)): the maximum local coordinate range
    """

    result = StgDomain.FieldVariable_GetMinAndMaxLocalCoords( fieldVar );
    return result[0], result[1]

def FieldVariable_GetMinAndMaxGlobalCoords( fieldVar ):
    """
    Returns the domain of this field variable globally as two tuples (min,max)
    
    Args:
        fieldVar (Swig FieldVariable*): Field Variable to query
    Returns:
        min (tuple(float)): the minimum global coordinate range
        max (tuple(float)): the maximum global coordinate range
    """

    result = StgDomain.FieldVariable_GetMinAndMaxGlobalCoords( fieldVar );
    return result[0], result[1]


def FeVariable_Integrate( feVar, gaussSwarm=None ):
    """
    Returns the integral of this FE variable.
    Note that this is currently only compatible with scalar fields.
    
    Args:
        fieldVar (Swig FeVariable*): Field Variable to query
        gaussSwarm (Swig Swarm*): (optional) Gauss swarm to integrate over. If nothing provided,
                                             will try and locate gauss swarm with name "gaussSwarm"
    Returns:
        integral (float): the integral of the FeVariable.
    """
    if gaussSwarm == None:
        gaussSwarm = GetLiveComponent("gaussSwarm")

    if feVar.fieldComponentCount != 1:
        print "Error: This function currently only supports scalar fields."
        return

    return StgFEM.FeVariable_Integrate( feVar, gaussSwarm )

def Swarm_GetVariables( swarm ):
    """
    Returns a python list of swarm variables. 
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        variables (list): python list of variables, variable data stored as tuples: (name, pointer, datatype)
    """
    varList = []
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    for i in range(0,variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes,0)
            # guyTup = (guy.name, guy, SwarmVariable_GetType(datatype))
            guyList = [guy.name, guy, SwarmVariable_GetType(datatype)]
            varList.append(guyList)
            #varList.append(guyTup)
    return varList



def Swarm_GetVariablesAsDict( swarm ):
    """
    Returns a python dictionary of swarm variables with variable names as keys. 
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        variables (dict): python dictionary of variables by variableName, variable data stored as list: (pointer, datatype)
    """

    varDict = {}
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    for i in range(0,variableCount):
        swarmObject = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if swarmObject.variable:
            datatype = StGermain.VariableTypeArrayDeref(swarmObject.variable.dataTypes,0)
            varDict[swarmObject.name] = {}
            varDict[swarmObject.name]['swarmVariable'] = swarmObject
            varDict[swarmObject.name]['dataType'] = SwarmVariable_GetType(datatype)
    
    return varDict

def Swarm_PrintVariables( swarm ):
    """
    Formatted print to standard out of swarm's variables.
    Note that some swarm variables (such as MaterialSwarmVariables) are excluded as the user cannot modify their data.
    
    Args:
        swarm (Swig Swarm*): Swarm
    Returns:
        Nothing
    """
    varNameList = []
    varTypeList = []
    maxLen = 0
    swarmreg = swarm.swarmVariable_Register
    variableCount = StGermain.Stg_ObjectList_CountFunc(swarmreg.objects)
    goodVarCount = 0
    for i in range(0,variableCount):
        guy = StGermain.Stg_ObjectList_AtFunc(swarmreg.objects, i)
        if guy.variable:
            goodVarCount = goodVarCount + 1
            varNameList.append(str(guy.name))
            if len(guy.name) > maxLen:
                maxLen = len(guy.name)
            datatype = StGermain.VariableTypeArrayDeref(guy.variable.dataTypes,0)
            varTypeList.append(SwarmVariable_GetType(datatype))
            #print "Name = %-40s Type = %-40s " % (guy.name, SwarmVariable_GetType(datatype))
    for i in range(0,goodVarCount):
        print "Name =", varNameList[i].ljust(maxLen+1), "Type =", varTypeList[i]

def SwarmVariable_GetType( datatype ):
    """
    Returns a string describing the variable type for the swarm variable datatype (as enum) provided
    
    Args:
        datatype (int / enum): swarm variable datatype enumeration
    Returns:
        variableType (str)
    """
    if   datatype == StGermain.Variable_DataType_Variable:
        return "Variable_DataType_Variable"
    elif datatype == StGermain.Variable_DataType_Char:
        return "Variable_DataType_Char"
    elif datatype == StGermain.Variable_DataType_Short:
        return "Variable_DataType_Short"
    elif datatype == StGermain.Variable_DataType_Int:
        return "Variable_DataType_Int"
    elif datatype == StGermain.Variable_DataType_Float:
        return "Variable_DataType_Float"
    elif datatype == StGermain.Variable_DataType_Double:
        return "Variable_DataType_Double"
    elif datatype == StGermain.Variable_DataType_Pointer:
        return "Variable_DataType_Pointer"
    elif datatype == StGermain.Variable_DataType_Size:
        return "Variable_DataType_Size"
    else:
        return "UNKNOWN"


def SwarmVariable_GetValueAt( swarmVar, localParticleIndex ):
    """
    Returns the swarm variable value from a local particle
    
    Args:
        swarmVar (Swig SwarmVariable*): Swarm variable to query.
        localParticleIndex (int)      : Local particle index for particle of interest
    Returns:
        result (list(swarmVar(type))) : The result as a list.  Value(s) in list are of the same type as the swarm variable.
    """
    toreturn = []
    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes,0)
    if   datatype == StGermain.Variable_DataType_Int:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtInt(    swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Float:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtFloat(  swarmVar.variable, localParticleIndex, ii ))
    elif datatype == StGermain.Variable_DataType_Double:
        for ii in range(0,swarmVar.dofCount):
            toreturn.append(StGermain.Variable_GetValueAtDouble( swarmVar.variable, localParticleIndex, ii ))
    else:
        return "Sorry, the swarm variable datatype", SwarmVariable_GetType(dataType), "is not supported."

    return toreturn


def SwarmVariable_SetValueAt( swarmVar, localParticleIndex, values ):
    """
    Sets the swarm variable value for a local particle
    
    Args:
        swarmVar (Swig SwarmVariable*): Swarm variable to set value for.
        localParticleIndex (int)      : Local particle index for particle of interest
        values (list)                 : values to set on particle
    Returns:
        None
    """
    import c_arrays
    
    if( swarmVar.dofCount < len(values) ):
        return "Error: size of values list is greater than variable dofCount"

    datatype = StGermain.VariableTypeArrayDeref(swarmVar.variable.dataTypes,0)

    if   datatype == StGermain.Variable_DataType_Int:
        valuePtr = c_arrays.IntArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Float:
        valuePtr = c_arrays.FloatArray(swarmVar.dofCount)
    elif datatype == StGermain.Variable_DataType_Double:
        valuePtr = c_arrays.DoubleArray(swarmVar.dofCount)
    else:
        return "Sorry, the swarm variable datatype", SwarmVariable_GetType(dataType), "is not supported."

    for ii in range(0,len(values)):
        valuePtr[ii] = values[ii]

    StGermain.Variable_SetValue( swarmVar.variable, localParticleIndex, valuePtr.cast() )


def NewComponentEntryInDictionary( name, Type, globalDict ):
    """
    Creates a new item in the components section of the dictionary

    Args:
        name (string)            : The name of the component (must be unique) 
        Type (string)             : The component Type
        globalDict (dictionary) : The runtime global dictionary

    Returns:
        dictionary node to be populated with the component data
    """

    if name == "":
        print "NewComponentEntryInDictionary requires a (unique) name"
        return

    if type == "":
        print "NewComponentEntryInDictionary requires a (valid) type"
        return
    
    newComponentDict = dict()
    globalDict["components"][name] = newComponentDict
    newComponentDict["name"] = name   # useful to have this ... 
    newComponentDict["Type"] = Type

    # and that's all we can do without knowing more about the component's type

    return newComponentDict


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



def StgComponentToGlobalDict( **cptDefinitionArgs ):

    """
    Updates the global python dictionary components with the information supplied. 

    Args:

        name=(string)            : The name of the component (must be unique) 
        Type=(string)             : The component type
        mergeType (string)        : The equivalent of the U/W mergeType (replace or merge)
        XXX=(some data)            : As defined by the component itself

    Returns:
        stgDictionary dictionary node for the (new) component 
    """

    if not "name" in cptDefinitionArgs:
        print "Error: NewComponentEntryInStgDict requires a (unique) 'name' "
        print "Cpt: ", cptDefinitionArgs
        return
    
    if not "Type" in cptDefinitionArgs:
        print "Error: NewComponentEntryInStgDict requires a (valid) 'Type' "
        print "Cpt: ", cptDefinitionArgs
        return

    globalDict = GetCurrentPythonDictionary()

    if "mergeType" in cptDefinitionArgs and cptDefinitionArgs["mergeType"] == "replace":
        del globalDict["components"][cptDefinitionArgs["name"]]
        globalDict["components"][cptDefinitionArgs["name"]] = cptDefinitionArgs
    else:
        if cptDefinitionArgs["name"] in globalDict["components"]:
            globalDict["components"][cptDefinitionArgs["name"]].update(cptDefinitionArgs)
        else:
            globalDict["components"][cptDefinitionArgs["name"]] = cptDefinitionArgs


    return globalDict["components"][cptDefinitionArgs["name"]]


# def StartVisualisationWebServer(host='127.0.0.1', port=-99999, certfile=None):
#     """
#     Starts a visualisation web server (via root processor) in a new thread.
#     An Init function first needs to be called so that MPI_Init has delegated processor Ids. 

#     Args:
#         host     (string)        : (optional) The server host address.  Defaults to localhost (127.0.0.1).
#         port     (   int)        : (optional) The server listening port.  Defaults to 8999 and autoincrementing if necessary.
#         certfile (string)       : (optional) A certificate file (including path) for running an HTTPS web server.

#     Returns:
#         Nothing
#     """
#     global _data

#     if _data == None:
#         print "StGermain has not been initialised yet. You need to run one of the Init functions."
#         return

#     if rank() == 0:
#         defaultport = 8999
#         # lets first halt any running instances
#         StopVisualisationWebServer()

#         if port < 0:
#             _port = defaultport
#         else:
#             _port = port

#         import SimpleHTTPServer
#         import SocketServer
#         import threading
#         global _httpd

#         SocketServer.TCPServer.allow_reuse_address = True
#         Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

#         if port < 0:
#             trying = True
#             while(trying):
#                 try:
#                     _httpd = SocketServer.TCPServer((host, _port), Handler)
#                     trying = False
#                 except Exception, e:
#                     _port +=1
#                     if _port > defaultport + 10:
#                         print "Error. Tried multiple ports without success. Perhaps try setting a port explicitly"
#                         print "using the port= option, or shutdown any already running servers."
#                         print str(e)
#                         raise
#         else:
#             try:
#                 _httpd = SocketServer.TCPServer((host, _port), Handler)
#             except Exception, e:
#                 print "Error. Unable to open port ", _port
#                 print str(e)
#                 raise

#         if certfile:
#             import ssl
#             _httpd.socket = ssl.wrap_socket(_httpd.socket, certfile=certfile, server_side=True)

#         httpd_thread = threading.Thread(target=_httpd.serve_forever)
#         httpd_thread.setDaemon(True)
#         httpd_thread.start()
#         print("Serving visualisation web server at %s:%i" % (host,_port))

# def StopVisualisationWebServer():
#     """
#     Halts any running visualisation servers

#     Args:
#         None

#     Returns:
#         Nothing
#     """
#     global _httpd
#     if _httpd:
#         _httpd.shutdown()
#         _httpd = None
#         print "Visualisation web server halted"


# ===========================    



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

    info = _globalDict["info"]
    
    return info

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
    NewComponentEntryInStgDict( _globalDict,
                                name = "context",
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

### Global Variables of the underworld module should go here. 
### If they are _variableName they remain private (but can be accessed by importing _underworld)
### If they are not _variableName they will be exposed through the import *, but it is difficult
### to use them as global variables due to issues of managing the scope. 


# init all 
_data = None
_isConstructed = False
_isBuilt = False
# _httpd = None

_globalDict = None

## Access global variables

def isConstructed():
    return _isConstructed

def isBuilt():
    return _isBuilt

def isInitialised():
    
    if (_data != None):
        return True
    else:
        return False

def data():
    return _data






