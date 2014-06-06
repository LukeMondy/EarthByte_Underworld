# Swarms package - configures swarms
import underworld as uw
##############################################################################
##  Adds Visualisation entries to the python dictionary  
##############################################################################

'''
  Adds Visualisation entries to the python dictionary 
'''
def titleCreate(name="randomTitle", string="Random Plot"):

    globalDict = uw.GetCurrentPythonDictionary()
    # Test if we already have a title of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    
    title = uw.NewComponentEntryInStgDict( globalDict, 
                                           name = name,
                                           Type = "lucTitle",
                                           string = str(string)
                                           )
    return title

def fieldVariableBorderCreate(name="", fieldVariable=""):
    """
    The fieldVariable can be an FeVariable or an OperatorFeVariable
    """
    globalDict = uw.GetCurrentPythonDictionary()
    if fieldVariable=="":
        uw.utils.sendError("Must provide an FeVariable or OperatorFeVariable name")
        return -1

    if name != "":
        # Test if we already have a fieldVariable border of the name we are going to use here.
        name = uw.utils.checkForNewComponentName(globalDict, name)
    else:
        name = fieldVariable+"Border"
    
    title = uw.NewComponentEntryInStgDict( globalDict, 
                                           name = name,
                                           Type = "lucFieldVariableBorder",
                                           FieldVariable = fieldVariable,
                                           )
    return title

    
def colourMapCreate(name="defaultColourMap",colours="Purple Blue Green Yellow Orange Red", dynamicRange=True, minimum=0.0, maximum=100.0):
    """
    Set colour map for plot.
    
    e.g.
    colours = "Purple Blue Green Yellow Orange Red"
      gives a colour map ranging from Purple through to Red matching highest to lowest values
      in the field being plotted
    """

    globalDict = uw.GetCurrentPythonDictionary()
    # Need to test if we already have a colour map of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)

    cmap = uw.NewComponentEntryInStgDict( globalDict, 
                                          name         = name,
                                          Type         = "lucColourMap", 
                                          colours      = colours,
                                          dynamicRange = str(dynamicRange),
                                          minimum = str(minimum),
                                          maximum = str(maximum)
                                          )
    return cmap

def colourBarCreate(name="defaultColourBar", colourMap=""):
    """
    Give a colour bar for the data range of a field plot using a colourMap.
    """

    globalDict = uw.GetCurrentPythonDictionary()
    # Need to test if we already have a colour bar of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    if colourMap=="":
        uw.utils.sendError("Must provide a ColourMap name")
        return -1

    cmap = uw.NewComponentEntryInStgDict( globalDict, 
                                          name      = name,
                                          Type      = "lucColourBar", 
                                          ColourMap = colourMap
                                          )
    return cmap


def scalarFieldMapCreate(name="", colourMap="", fieldVariable="", resolution=128):
    """
    Create a scalar field plot.

    The fieldVariable can be an FeVariable or an OperatorFeVariable
    """
    globalDict = uw.GetCurrentPythonDictionary()

    if fieldVariable=="":
        uw.utils.sendError("Must provide an FeVariable or OperatorFeVariable name")
        return -1

    if colourMap=="":
        # just create a default map
        colourMapDict = colourMapCreate(name=fieldVariable+"ColourMap",colours="Purple Blue Green Yellow Orange Red", dynamicRange="True")
        colourMap = colourMapDict["name"]
    if name != "":
        # Test if we already have a colour map of the name we are going to use here.
        name = uw.utils.checkForNewComponentName(globalDict, name)
    else:
        name = fieldVariable+"ScalarFieldMap"

    sfield = uw.NewComponentEntryInStgDict( globalDict, 
                                            name          = name,
                                            Type          = "lucScalarField", 
                                            FieldVariable = fieldVariable,
                                            ColourMap     = colourMap,
                                            resolution    =str(resolution)
                                            )
    return sfield

def scalarContoursCreate(name="", colour="black", fieldVariable="", interval=0.2, resolution=20, lineWidth=2, showValues="False"):
    """
    Create a scalar field plot.

    The fieldVariable can be an FeVariable or an OperatorFeVariable
    """
    globalDict = uw.GetCurrentPythonDictionary()

    if fieldVariable=="":
        uw.utils.sendError("Must provide an FeVariable or OperatorFeVariable name")
        return -1

    if name != "":
        # Test if we already have a name we are going to use here.
        name = uw.utils.checkForNewComponentName(globalDict, name)
    else:
        name = fieldVariable+"Contours"

    sfield = uw.NewComponentEntryInStgDict( globalDict, 
                                            name             = name,
                                            Type              = "lucScalarField", 
                                            FieldVariable = fieldVariable,
                                            resolution      = str(resolution),
                                            colour            = colour,
                                            interval          =  str(interval),
                                            lineWidth       = str(lineWidth),
                                            showValues  = str(showValues)
                                            )
    return sfield

def swarmViewerCreate(name="", colourMap="", pointSize=2.0, swarm="", variable="colour"):
    """
    Create a swarm viewer.
    """
    globalDict = uw.GetCurrentPythonDictionary()

    if swarm=="":
        uw.utils.sendError("Must provide a Swarm name")
        return -1

    if colourMap=="":
        # look for a default colour map; create if doesn't exist
        if "defaultColourMap" in globalDict["components"]:
            colourMap="defaultColourMap"
        else:
            colourMapDict = colourMapCreate(name="defaultColourMap",colours="Purple Blue Green Yellow Orange Red", dynamicRange="True")
            colourMap = colourMapDict["name"]
    if name != "":
        # Test if we already have the name we are going to use here.
        name = uw.utils.checkForNewComponentName(globalDict, name)
    else:
        name = swarm+"Viewer"

    sviewer = uw.NewComponentEntryInStgDict( globalDict, 
                                             name       = name,
                                             Type       = "lucSwarmViewer", 
                                             Swarm      = swarm,
                                             ColourVariable = swarm+"-"+variable,  # I think this guy gets created in code so is not in dictionary
                                             ColourMap  = colourMap,
                                             pointSize  = str(pointSize)
                                             )
    return sviewer

def cameraCreate(name="camera", coordZ=1.5, centreFieldVariable=""):
    """
    Create a camera

    """
    globalDict = uw.GetCurrentPythonDictionary()

    if centreFieldVariable=="":
        uw.utils.sendError("Must provide an FeVariable or OperatorFeVariable name")
        return -1

    # Test if we already have a camera  of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)

    sfield = uw.NewComponentEntryInStgDict( globalDict, 
                                            name          = name,
                                            Type          = "lucCamera", 
                                            CentreFieldVariable = centreFieldVariable,
                                            coordZ        = str(coordZ)
                                            )
    return sfield

# camera lucColourBar lucScalarField lucTitle
def viewPortCreate(name="", titleEntry="", colourBar="", scalarFieldMap="", fieldVariable="", camera="", border="", compositeEachObject="false", **extra):
    """
    Create a view port

    """
    # still have to make this colourMap border stuff a little clearer
    globalDict = uw.GetCurrentPythonDictionary()
    if fieldVariable=="" and scalarFieldMap=="":
        uw.utils.sendError("Must provide either an FeVariable/OperatorFeVariable name for fieldVariable\n"+ 
                           "or a scalarFieldMap map for scalarFieldMap")
        return -1
    if scalarFieldMap=="": # Then fieldVariable is not blank
        # create scalarFieldMap which also needs a colourMap
        # scalarFieldMapCreate can make its own defaultColourMap but let's create it and name it based on the fieldVariable
        cmap = colourMapCreate( name = fieldVariable+"ColourMap")
        sfield = scalarFieldMapCreate( fieldVariable = fieldVariable, colourMap = cmap["name"])
        scalarFieldMap=sfield["name"]
        # if we have made a new colourMap then any colourBar passed in doesn't make sense?
        cbar = colourBarCreate( colourMap = cmap["name"], name = fieldVariable+"ColourBar" )
        colourBar = cbar["name"]
        if titleEntry=="":
            title = titleCreate(name=fieldVariable+"Title", string = fieldVariable)
            titleEntry = title["name"]
        if name=="":
            name = fieldVariable+"VP"
        if border=="":
            bord = borderCreate(name="border", fieldVariable=fieldVariable)
            border=bord["name"]
    if fieldVariable=="": # then scalarFieldMap is not blank
        # if no colourBar then let's not worry about it?
        if titleEntry=="":
            title = titleCreate(name=scalarFieldMap+"Title", string = scalarFieldMap)
            titleEntry = title["name"]
        if name=="":
            name = scalarFieldMap+"VP"

    if camera=="":
        # just create a camera using the fieldVariable?
        cam=cameraCreate(name="camera", coordZ=1.5, centreFieldVariable=fieldVariable)
        camera=cam["name"]


    drawingObjectList=[]
    drawingObjectList.append(titleEntry)
    if colourBar!="":
        drawingObjectList.append(colourBar)
    drawingObjectList.append(scalarFieldMap)
    drawingObjectList.append(border)
    
    for key in extra.keys():
        drawingObjectList.append(extra[key])

    vport = uw.NewComponentEntryInStgDict( globalDict, 
                                           name          = name,
                                           Type          = "lucViewport", 
                                           Camera        = camera,
                                           compositeEachObject = str(compositeEachObject),
                                           DrawingObject = drawingObjectList
                                           )

    return vport



def databaseCreate(name="database", singleFile="False", 
                   writeEvery=1, splitTransactions="True",
                   filename="gLuciferDatabase",
                   drawingObjectList=[]):
    
    globalDict = uw.GetCurrentPythonDictionary()

    # Need to test if we already have a database of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    
    dbObject = uw.NewComponentEntryInStgDict( globalDict, 
                                              name      = name,
                                              Type      = "lucDatabase",
                                              singleFile=singleFile,
                                              writeEvery=writeEvery,
                                              filename=filename,
                                              DrawingObject = drawingObjectList
                                              )
    return dbObject

def windowCreate(name="window", viewPortList="", database="database"):
    
    globalDict = uw.GetCurrentPythonDictionary()
    # Need to test if we already have a database of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    
    if viewPortList=="":
        uw.utils.sendError("Must provide viewPorts as list of strings")
        return -1

    window = uw.NewComponentEntryInStgDict( globalDict, 
                                            name     = name,
                                            Type     = "lucWindow",
                                            Viewport = viewPortList,
                                            Database = database
                                            )
    
    return window

def borderCreate(name="border", fieldVariable=""):
 
    globalDict = uw.GetCurrentPythonDictionary()
    # Need to test if we already have a database of the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    
    if name=="":
        name = fieldVariable+"Border"

    border = uw.NewComponentEntryInStgDict( globalDict, 
                                            name     = name,
                                            Type     = "lucFieldVariableBorder",
                                            FieldVariable = fieldVariable
                                            )

    return border


def fieldArrowsCreate(name="fieldArrows", fieldVariable="", arrowHeadSize=0.15, lengthScale=0.15, colour="black"):
    """
    Create vector field Arrows
    """
    globalDict = uw.GetCurrentPythonDictionary()
    # Test if we already have the name we are going to use here.
    name = uw.utils.checkForNewComponentName(globalDict, name)
    
    if fieldVariable=="":
        uw.utils.sendError("Must provide an FeVariable or OperatorFeVariable name")
        return -1

    if name=="":
        name = fieldVariable+"Arrows"
    
    arrows = uw.NewComponentEntryInStgDict( globalDict, 
                                            name          = name,
                                            Type          = "lucVectorArrows", 
                                            VectorVariable = fieldVariable,
                                            Colour        = colour,
                                            arrowHeadSize = str(arrowHeadSize),
                                            lengthScale   = str(lengthScale)
                                            )
    return arrows
