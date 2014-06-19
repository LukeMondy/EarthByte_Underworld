import sys as _sys
import _stgermain
import libUnderworld
import dictionary

### These functions allow us to build, merge and manipulate the global dictionary as a python object
### and, later, it can be converted back to XML strings for the blessed StGermain to injest.     


def Init(args=[], prependCommandLineArgs=True):
    """
    Runs the Underworld init phase using the arguments provided at the script command line.

    Args:
        args (list, str): List of extra arguments, provided as a list or as a string.  Default is empty.
        prependCommandLineArgs (bool):  Should we prepend any arguments passed in at the command line?  Default is True.
    Returns:
        Nothing
    """

    Finalise()

    if type(args) == str:
        args = args.split()

    if prependCommandLineArgs:
        args = _sys.argv[0:]+args 
    else:
        args = _sys.argv[0:1]+args

    _stgermain.StgInit( args )

    uwdict = _stgermain.GetStgDictionaryAsPyDict(_stgermain.getData().dictionary)

    uwdict = dictionary.GetConformingDictionary(uwdict)

    if "parameters" not in uwdict.keys():
        uwdict["parameters"] = {}
    if "info" not in uwdict.keys():
        uwdict["info"] = {}

    dictionary.SetDictionary( uwdict )
      
    return


def Construct( buildNow=True ):
    """
    Runs the Underworld Construct simulation phase.
    An Init function must be called before this.

    Args:
        buildNow (bool):  Should the build stage be executed immediately. Default is True.
    Returns:
        Nothing
    """
    _AssertConstructStateIs(False)

    globalDict = dictionary.GetDictionary()
    
    _stgermain.StgConstruct(globalDict, setAsRootDict=True)
    _setConstructFlag(True)
    _stgermain.WriteFlattenedFile(globalDict)
    
    # lets go right ahead and build
    if buildNow:
       _stgermain.StgBuild(globalDict)
       _stgermain.StgInitialise(globalDict)

def RunMainLoop(context=None):
    """
    Runs the Underworld main simulation loop.
    An Init function must be called before this, as well as the Construct phase.

    Args:
        context (Swig Context*): Context for current simulation.  Default is none.
    Returns:
        Nothing
    """

    _AssertConstructStateIs(True)

    if context == None:
        context = _stgermain.GetLiveComponent("context")
    assert context, "No context found. You must either pass in a context, or one must exists in the dictionary."

    libUnderworld.StGermain.Stg_Component_Execute( context, None, True)

def Step(context=None, steps=1):
    """
    Performs required number of steps of the simulation. 
    An Init function must be called before this, as well as the Construc routine.
    This step will ignore all simulation flow control parameters (such as maxTimeStep).

    Args:
        context (Swig Context*): (optional) Context for current simulation.  
            If none provided, will search LiveComponentRegister    for component with name "context".
        steps (int): is the number of steps to take
    Returns:
        Nothing
    """
    _AssertConstructStateIs(True)

    if context == None:
       context = _stgermain.GetLiveComponent("context")
    assert context, "No context found. You must either pass in a context, or one must exists in the dictionary."

    for i in range(0,steps):
        if( context.needUpdate == 1 ):
            libUnderworld.StGermain.AbstractContext_Update(context)

        libUnderworld.StGermain.AbstractContext_Step(context, context.dt)

def Finalise():
    """
      Finalise the underworld simulation and tools.  Clears the global dictionary.
      
      Args:
         None
    """
    _stgermain.StgFinalise()
    _resetAll()

def _resetAll():
    global _isConstructed
    _setConstructFlag( False )
    dictionary.ClearDictionary()


def _AssertConstructStateIs(isConstructed):
   if isConstructed:
      assert ConstructState() is isConstructed, "Underworld has not been constructed.  You will need to call Construct() before you can perform this operation."
   else:
      assert ConstructState() is isConstructed, "Underworld has already been constructed.  You will need to call Init() again before you can perform this operation."

def ConstructState():
   """
      Returns the construct state of the simulation.  True=Constructed.
      
      Args:
      None
   """
   global _constructFlag
   return _constructFlag

def _setConstructFlag(isConstructed):
   global _constructFlag
   _constructFlag = isConstructed


# init all 
_constructFlag = False

