import sys as _sys
import os as _os
import tempfile as _tempfile
import subprocess as _subprocess
import xml.etree.cElementTree as _ET
import collections as _collections
from libUnderworld import *
import utils
import weakref
import abc

class StgStaleComponent(object):
    def __getattr__(self, name):
        raise AttributeError("""Error. This object is now stale and cannot be used. All existing objects of (sub)type 'StgCompoundComponent' become stale when Finalise() or Init() are called.""")

class StgCompoundComponent(object):
    """ 
    This class ties multiple StGermain components together into a single python object.
    The life cycle of the objects (construction/build/destruction) are handled automatically  
    
    """
    __metaclass__ = abc.ABCMeta
    
    _livingInstances = []
    
    def __new__(cls, objectDict, *args, **kwargs):
        """
        New function. Creates stgermain instances of underlying objects.
        
        Args:
            objectDict (dict)   : Objects dictionary to be provided by child.  Specifies object name (key) and object type (value).

        Returns:
            New created instance of child class.
        
        """
        if not isinstance(objectDict, dict):
            raise TypeError("object passed in must be of python type 'dict' or subclass")

        # lets go ahead and create the python instance of this object
        self = object.__new__(cls, *args, **kwargs)

        import string
        import random
        self._id = "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

        # ok, create stgermain objects
        # first rearrange to create stg compatible dictionary
        newObjDict = {}
        for compName, compType in objectDict.iteritems():
            newObjDict[self._id + "_" + compName] = { "Type": compType }
        fullDictionary = {"components": newObjDict}
        # create
        pointerDict = StgCreateInstances(fullDictionary)
        for objName, objPointer in pointerDict.iteritems():
            # note here we strip out the id part to make it user friendly
            setattr(self, objName.replace(self._id,"_clib"), objPointer)
                
        return self

    def __init__(self, modulesList = None):
        """ 
        Initialisation function.  All components in provided dictionary are constructed, build & initialised here.
        
        Args:
            modulesList (list)   : List of StGermain modules to load. Default = None
        """
        if modulesList:
            if not isinstance(modulesList, list):
                raise TypeError("object passed in must be of python type 'list' or subclass")
            LoadModules( {"import":modulesList} )
    
        self.componentDictionary = _collections.OrderedDict()
        # ok... let child classes fill dictionary.
        self._addToStgDict()
        self.fullDictionary = {"components": self.componentDictionary}
        StgConstruct(self.fullDictionary)
        StgBuild(self.fullDictionary)
        StgInitialise(self.fullDictionary)
        self._isAlive = True
        self._weakref = weakref.ref(self)
        StgCompoundComponent._livingInstances.append(self._weakref)
        
        super(StgCompoundComponent, self).__init__()

        # use the following dictionary to store user friendly names of stg components in child classes

    def __del__(self):
        """ 
        Destructor method.  Calls destroy & delete phases.
        
        """
        # check if we have the '_isAlive' attribute, incase one of the child classes died during __init__
        if hasattr(self,"_isAlive") and self._isAlive:
            StgDestroy(self.componentDictionary)
            StgDelete(self.componentDictionary)
            StgCompoundComponent._livingInstances.remove(self._weakref)
            self._isAlive = False

    @abc.abstractmethod
    def _addToStgDict(self):
        """ This function needs to be set by child class.
            It allows each child class to enter values into the component dictionary, and then call the parent _addToStgDict method.
            
            Args:
                None
            Returns:
                None
        """

    def _setterAssertNotConstructed(self):
        if hasattr(self,"_isAlive"):
            raise RuntimeError("This object attribute can only be set via the object constructor.")




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
                if dicttag in ['struct', 'StGermainData']:
                    if not childkey in value:
                        # if no entry with this key, add one
                        value[childkey] = valueGuy
                elif dicttag == 'list':
                    # if entry already exists,
                    try:
                        # assume it's a list, and add item to list
                        value[childkey].append(valueGuy)
                    except:
                        # if that fails, lets create a list and add item
                        value[childkey] = []
                        value[childkey].append(valueGuy)
                else:
                    raise TypeError("Param cannot contain sub values. key={}, value={}".format(childkey,valueguy))
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
    elif string.lower() == "true":
        return True
    elif string.lower() == "false":
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
    for k, v in inputDict.items():
        _itemToElement(v, k, root)

    return root


def _itemToElement(inputItem, inputItemName, inputEl):
    if type(inputItem) == list:
        subEl = _ET.SubElement(inputEl, 'list')
        if inputItemName != '':
            subEl.attrib['name'] = inputItemName
        for item in inputItem:
            _itemToElement(item, '', subEl)
    elif type(inputItem) in [_collections.OrderedDict, dict]:
        subEl = _ET.SubElement(inputEl, 'struct')
        if inputItemName != '':
            subEl.attrib['name'] = inputItemName
        for k, v in inputItem.items():
            _itemToElement(v, k, subEl)
    elif type(inputItem) in [str, float, int, bool, unicode]:
        subEl = _ET.SubElement(inputEl, 'param')
        if inputItemName != '':
            subEl.attrib['name'] = inputItemName
        subEl.text = str(inputItem)
    elif not inputItem:
        subEl = _ET.SubElement(inputEl, 'param')
        if inputItemName != '':
            subEl.attrib['name'] = inputItemName
        subEl.text = "\t"
    else:
        raise TypeError("Unknown type encountered. key={}, value={}".format(inputItemName,inputItem))


def _indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def GetStgDictionaryAsPyDict( stgDict ):
    """
    Gets a copy of the current StGermain dictionary
    This routines can only be utilised after the Init phase has been completed.

    Args:
        stgDict (Swig StGermain Dictionary*) : input dictionary to convert to python dictionary.
    Returns:
        dict (OrderedDict): Python Ordered dictionary of Underworld components.

    """
    _AssertInitStateIs(True)
    ioHandler = StGermain.XML_IO_Handler_New()
    xmlString = StGermain._XML_IO_Handler_WriteAllMem(ioHandler, stgDict, None )
    StGermain.Stg_Class_Delete( ioHandler )
    root = _ET.fromstring(xmlString)
    dict = _elementToDict(root)

    return dict


def SetStgDictionaryFromPyDict( pyDict, stgDict ):
    """
       Sets the provided python dictionary as the StGermain dictionary.
       This routines can only be utilised after the Init phase has been completed, and
       should usually be run before the Construct phase.

       Args:
       pyDict (dict): Python dictionary to build a StGermain dictionary from.
       stgDict (Swig StGermain Dictionary*):  Pointer to StGermain dictionary to add python dictionary contents to.

       Returns:
       Nothing.
       """
    _AssertInitStateIs(True)
    root = _dictToUWElementTree(pyDict)
    xmlString = _ET.tostring(root, encoding = 'utf-8', method = 'xml')
    ioHandler = StGermain.XML_IO_Handler_New()
    StGermain.IO_Handler_ReadAllFromBuffer( ioHandler, xmlString, stgDict, "" )
    StGermain.Stg_Class_Delete( ioHandler )

    return


def LoadModules( pyUWDict ):
    """
       Loads any Toolboxes found within provided dictionary.

       Args:
           pyUWDict (dict): Python version of underworld root dictionary.

       Returns:
           Nothing.
    """
    _AssertInitStateIs(True)
    stgRootDict = StGermain.Dictionary_New()
    SetStgDictionaryFromPyDict( pyUWDict, stgRootDict )
    StGermain.ModulesManager_Load( StGermain.GetToolboxManagerInstance(), stgRootDict, "" )
    StGermain.Stg_Class_Delete( stgRootDict )


def StgInit( args=[] ):
    """
       Calls the StGermain Init function with provided arguments.

       Args:
           args (list): List of arguments to pass to Init stage.  Usually command line like arguments.  Default none.

       Returns:
           Nothing.
    """
    
    if getData():
        StgFinalise()
    setData( StGermain_Tools.StgInit( args ) )


def WriteFlattenedFile( uwdict, timestamp=None ):
    """
       Writes a flattened XML file using provided dictionary.  File generated is input.xml.

       Args:
       uwdict (dict): Dictionary to write to flattened xml file.
       timestamp (str):  If provided, generates timestampted filenames alongside input.xml.

       Returns:
       Nothing.
    """
    _AssertInitStateIs(True)
    if not isinstance(uwdict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")
    stgDict = StGermain.Dictionary_New()
    SetStgDictionaryFromPyDict( uwdict, stgDict )
    StGermain.stgGenerateFlattenedXML( stgDict, None, None )
    StGermain.Stg_Class_Delete(stgDict)

def StgCreateInstances( pyUWDict ):
    """
        Creates instances for all components within pyUWDict.
        
        Args:
        pyUWDict (dict): Underworld root type dictionary containing components and plugins.
        
        Returns:
        pointerDict (dict): Dictionary mapping component names to stg pointer.
    """
    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    LoadModules(pyUWDict)

    stgRootDict = StGermain.Dictionary_New()
    SetStgDictionaryFromPyDict( pyUWDict, stgRootDict )

    stgCompDict = StGermain.Dictionary_Entry_Value_AsDictionary( StGermain.Dictionary_Get( stgRootDict, "components" ) )

    cf = StGermain.Stg_ComponentFactory_New( stgRootDict, stgCompDict )

    # lets create instances of components
    StGermain.Stg_ComponentFactory_CreateComponents( cf )

    StGermain.Stg_Class_Delete(cf)
    StGermain.Stg_Class_Delete(stgRootDict)

    pointerDict = {}
    # lets go ahead and construct component
    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), compName )
            pointerDict[compName] = compPointer
# disable this as its broken because i believe the context creates the plugins during assignfromxml
#    if "plugins" in pyUWDict:
#        for guy in pyUWDict["plugins"]:
#            print guy["Type"]
#            compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), guy["Type"] )
#            pointerDict[guy["Type"]] = compPointer

    return pointerDict


def StgConstruct( pyUWDict, setAsRootDict=False):
    """
       Calls the construct phase for all components & plugins found in provided dictionary.

       Args:
         pyUWDict (dict): Underworld root type dictionary containing components and plugins.
         setAsRootDict   (bool):  If true, retains generated stg_componentfactory and stg root dictionary. Default false.

       Returns:
         Nothing
    """
    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    stgRootDict = StGermain.Dictionary_New()
    SetStgDictionaryFromPyDict( pyUWDict, stgRootDict )
    # now lets de-alias
    StGermain.DictionaryUtils_AliasDereferenceDictionary( stgRootDict )

    stgCompDict = StGermain.Dictionary_Entry_Value_AsDictionary( StGermain.Dictionary_Get( stgRootDict, "components" ) )

    cf = StGermain.Stg_ComponentFactory_New( stgRootDict, stgCompDict )

    #
    # It is currently necessary to ensure that the context gets built/constructed first
    # It is currently possible to have components trying to build themselves with unitialised data
    # from other components if the context has not been built already.
    #
    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            try:
                compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), compName )
                if "context" in compName:
                    StGermain.Stg_Component_AssignFromXML( compPointer, cf, None, False )
            except:
                utils.sendWarning("Component \'%s\' not found in the live component register." % compName)
    
    
    # lets go ahead and construct component
    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            try:
                compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), compName )
                StGermain.Stg_Component_AssignFromXML( compPointer, cf, None, False )
            except:
                utils.sendWarning("Component \'%s\' not found in the live component register." % compName)

    if "plugins" in pyUWDict:
        for guy in pyUWDict["plugins"]:
            try:
                compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), guy["Type"] )
                #print(compPointer)
                StGermain.Stg_Component_AssignFromXML( compPointer, cf, None, False )
            except:
                utils.sendWarning("Component \'%s\' not found in the live component register in Construct phase!" % guy["Type"])

    # don't like this, but not much choice at this point
    # we retain the concept of a root dict, and also the component factor object, as some things require it during the build phase (annoyingly)
    # note that below is mainly book keeping / mem management.  the items (except compdict) are stored on the context in abstractcontext_assignfromxml()
    if setAsRootDict is True:
        StGermain.Stg_Class_Delete(getData().dictionary)
        getData().dictionary = stgRootDict
        StGermain.Stg_Class_Delete(getData().cf)
        getData().cf = cf
    else:
        StGermain.Stg_Class_Delete(cf)
        StGermain.Stg_Class_Delete(stgRootDict)


def StgBuild( pyUWDict=None ):
    """
       Calls the build phase for all components & plugins found in provided dictionary.

       Args:
       pyUWDict (dict): Underworld root type dictionary containing components and plugins.

       Returns:
       Nothing.
    """

    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            try:
                compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), compName )
                StGermain.Stg_Component_Build( compPointer, None, False )
            except:
                utils.sendWarning("Component \'%s\' not found in the live component register." % compName)
    if "plugins" in pyUWDict:
        for guy in pyUWDict["plugins"]:
            try:
                compPointer = StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), guy["Type"] )
                StGermain.Stg_Component_Build( compPointer, None, False )
            except:
                utils.sendWarning("Component \'%s\' not found in the live component register in Build phase!" % guy["Type"])

    # now run the original c-code Build,
    # by looping over the LiveComponentRegister and Build each component
    lcr = StGermain.LiveComponentRegister_GetLiveComponentRegister()
    for ii in range(0,StGermain.LiveComponentRegister_GetCount(lcr)):
        component = StGermain.LiveComponentRegister_At( lcr, ii );
        StGermain.Stg_Component_Build( component, None, False )

def StgInitialise( pyUWDict ):
    """
       Calls the Initialise phase for all components & plugins found in provided dictionary.

       Args:
       pyUWDict (dict): Underworld root type dictionary containing components and plugins.

       Returns:
       Nothing.
    """

    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            compPointer = GetLiveComponent( compName )
            StGermain.Stg_Component_Initialise( compPointer, None, False )
    if "plugins" in pyUWDict:
        for guy in pyUWDict["plugins"]:
            compPointer = GetLiveComponent( guy["Type"] )
            StGermain.Stg_Component_Initialise( compPointer, None, False )

    # now run the original c-code Initialise,
    # by looping over the LiveComponentRegister and Initialise each component
    lcr = StGermain.LiveComponentRegister_GetLiveComponentRegister()
    for ii in range(0,StGermain.LiveComponentRegister_GetCount(lcr)):
        component = StGermain.LiveComponentRegister_At( lcr, ii );
        StGermain.Stg_Component_Build( component, None, False )


def StgDestroy( pyUWDict ):
    """
       Calls the Destroy phase for all components & plugins found in provided dictionary.

       Args:
       pyUWDict (dict): Underworld root type dictionary containing components and plugins.

       Returns:
       Nothing.
    """
    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")

    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            compPointer = GetLiveComponent( compName )
            StGermain.Stg_Component_Destroy( compPointer, None, False )
    if "plugins" in pyUWDict:
        for guy in pyUWDict["plugins"]:
            compPointer = GetLiveComponent( guy["Type"] )
            StGermain.Stg_Component_Destroy( compPointer, None, False )

def StgDelete( pyUWDict ):
    """
        Calls the Delete phase for all components & plugins found in provided dictionary.
        
        Args:
        pyUWDict (dict): Underworld root type dictionary containing components and plugins.
        
        Returns:
        Nothing.
        """
    _AssertInitStateIs(True)
    if not isinstance(pyUWDict, dict):
        raise TypeError("object passed in must be of python type 'dict' or subclass")
    
    if "components" in pyUWDict:
        for compName, compDict in pyUWDict["components"].iteritems():
            compPointer = GetLiveComponent( compName )
            StGermain.Stg_Class_Delete( compPointer )
    if "plugins" in pyUWDict:
        for guy in pyUWDict["plugins"]:
            compPointer = GetLiveComponent( guy["Type"] )
            StGermain.Stg_Class_Delete( compPointer )


def StgFinalise():
    """
       Finalises / tears down the StGermain simulation.

       Args:
       None
       Returns:
       Nothing
    """
    if getData():
        # first delete any python based living components, as they will become invalide once finalise runs.
        # we iterate backwards as this makes most sense (delete objects in opposite order of creation.
        # note that the python objects themselves will remain (until they go out of scope), but the underlying
        # stgermain data will be destroyed & deleted
        while len(StgCompoundComponent._livingInstances) is not 0:
            lastGuyPos = len(StgCompoundComponent._livingInstances) - 1
            lastGuy = StgCompoundComponent._livingInstances[lastGuyPos]()
            lastGuy.__del__()
            lastGuy.__class__ = StgStaleComponent
        StGermain_Tools.StgFinalise( getData() )
        setData(None)
    return


def GetLiveComponent(compName):
    """
       Returns component with provided name if found within live component register.  Otherwise returns None.

       Args:
       compName (str):  Name of component to return.
       Returns:
       component (Swig Ptr):  Returns a pointer to the component object.  If not found, returns None.
    """
    _AssertInitStateIs(True)
    if not isinstance(compName, str):
        raise TypeError("object passed in must be of python type 'str' or subclass")

    try:
        return StGermain.LiveComponentRegister_Get( StGermain.LiveComponentRegister_GetLiveComponentRegister(), compName )
    except:
        print "Component \'%s\' not found in the live component register." % compName
        return None


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
            flattenXML = _os.path.join(build_path, "bin", "FlattenXML")
            current_dir = _os.getcwd()
            tempdir = _tempfile.mkdtemp()
            _os.chdir(tempdir)
            _subprocess.call([flattenXML, xmlFile])
            xmlFile = _os.path.join(tempdir, "output.xml")
            _os.chdir(current_dir)

    if type(xmlFile) == str:
        theFile = open(xmlFile, 'r')
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

    # remove any temp files ...
    # HERE

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
        theFile = open(jsonFile, 'w')
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
        theFile = open(jsonFile, 'r')
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


def _AssertInitStateIs(isInit):
    if isInit:
        assert getData(), "StGermain has not been initialised.  You will need to call StgInit() before you can perform this operation."
    else:
        assert not getData(), "StGermain has been initialised.  You will need to call StgFinalise() before you can perform this operation."


def getData():
    global _data
    return _data


def setData(data):
    global _data
    _data = data
    if _data:
        _data.cf = None

_data = None
