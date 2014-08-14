import underworld
import errno
import weakref
import underworld._stgermain as _stgermain
import _drawing

# lets create somewhere to dump data for this session
import os
try:
    tmpdir = os.environ['TMPDIR']
except:
    tmpdir = "/tmp"
os.path.join(tmpdir,"glucifer")
try:
    os.makedirs(tmpdir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

class Figure(_stgermain.StgCompoundComponent):
    """  The Figure class provides the base container for gLucifer drawing objects, and associated routines for image generation.
         Generally the user will use the figure() routine to generate or recall Figure instances.
         
         A minimal workflow might be as follows:
            fig = plt.figure()                  # generate a Figure instance
            fig.Surface(field="PressureField")  # add Surface drawing object
            fig.show()                          # show the image (in ipython notebook)
            
         Currently, the Surface, Points & VectorArrows are the supported drawing objects.  See help() on these objects for their respective options.
    """
    def __init__(self, num=None, figsize=(640,480), facecolour="white", edgecolour="white", title=None, axis=False, **kwargs):
        """ The initialiser takes as arguments 'num', 'figsize', 'facecolour', 'edgecolour', 'title' and 'axis'.   See help(Figure) for full details on these options.
        """
        if num and not isinstance(num,(str,int)):
            raise TypeError("'num' object passed in must be of python type 'str' or 'int'")
        if num and isinstance(num,str) and (" " in num):
            raise ValueError("'num' object passed in must not contain any spaces.")
        self._num = num

        if not isinstance(figsize,tuple):
            raise TypeError("'figsize' object passed in must be of python type 'tuple'")
        self._figsize = figsize

        if not isinstance(facecolour,str):
            raise TypeError("'facecolour' object passed in must be of python type 'str'")
        self._facecolour = facecolour

        if not isinstance(edgecolour,str):
            raise TypeError("'edgecolour' object passed in must be of python type 'str'")
        self._edgecolour = edgecolour

        if title and not isinstance(title,str):
            raise TypeError("'title' object passed in must be of python type 'str'")
        self._title = title

        if not isinstance(axis,bool):
            raise TypeError("'axis' object passed in must be of python type 'bool'")
        self._axis = axis

        self._drawingObjects = []
        
        super(Figure,self).__init__(**kwargs)


    def _addToStgDict(self):
        # lets build up component dictionary
        # append random string to provided name to ensure unique component names
        # call parents method
        super(Figure,self)._addToStgDict()

        uniqueid = self._getUniqueName(prefix=self.num)
        self._localNames[ "db"] = "db_"  + uniqueid
        self._localNames["win"] = "win_" + uniqueid
        self._localNames[ "vp"] = "vp_"  + uniqueid
        self._localNames["cam"] = "cam_" + uniqueid


        self.componentDictionary[self._localNames["db"]] = {
                            "Type"              :"lucDatabase",
                            "filename"          :"gluciferDB"+uniqueid,
                            "blocking"          :True,
                            "dbPath"            :tmpdir
        }
        self.componentDictionary[self._localNames["win"]] = {
                            "Type"              :"lucWindow",
                            "Database"          :self._localNames[ "db"],
                            "Viewport"          :[self._localNames["vp"]],
                            "width"             :self.figsize[0],
                            "height"            :self.figsize[1],
                            "backgroundColour"  :self.facecolour,
                            "useModelBounds"    :False
        }
        self.componentDictionary[self._localNames["vp"]] = {
                            "Type"              :"lucViewport",
                            "Camera"            :self._localNames["cam"],
                            "borderColour"      :self.edgecolour,
                            "border"            :1,
                            "title"             :self.title,
                            "axis"              :self.axis
        }
        self.componentDictionary[self._localNames["cam"]] = {
                            "Type"              :"lucCamera",
                            "useBoundingBox"    : True
        }


    @property
    def num(self):
        """    num (str,int): integer or string figure identifier. Must not contain spaces. optional, default: none
        """
        return self._num

    @property
    def figsize(self):
        """    figsize (tuple(int,int)): size of window in pixels, default: (640,480)
        """
        return self._figsize

    @property
    def facecolour(self):
        """    facecolour : colour of face background, default: white
        """
        return self._facecolour

    @property
    def edgecolour(self):
        """    edgecolour : colour of figure border, default: white
        """
        return self._edgecolour

    @property
    def title(self):
        """    title : a title for the image, default: None
        """
        return self._title

    @property
    def axis(self):
        """    axis : Axis enabled if true.  Default False.
        """
        return self._axis

    @property
    def drawingObjects(self):
        """    drawingObjects : list of objects to be drawn within the figure.
        """
        return self._drawingObjects

    def show(self):
        """    Shows the generated image inline within an ipython notebook
        
                Args:
                    None
                    
                Returns:
                    Ipython Image object (will be displayed inline in an ipython notebook).
                    If Ipython is not found, nothing is returned.
                    
        """
        
        try:
            from IPython.display import Image
            self._generateDB()
            self._generateImage()
            return Image(filename=self._findGeneratedFile())
        raise:
            pass

    def _findGeneratedFile(self):
        # lets determine what we are outputting (jpg,png)
        foundFile = None
        for extension in ("jpg", "jpeg", "png"):
            fname = os.path.join(tmpdir,self._localNames["win"]+".00000."+extension)
            if os.path.isfile(fname):
                foundFile = fname
                break
        
        if not foundFile:
            raise RuntimeError("The required rendered image does not appear to have been created. Please contact developers.")
        
        return os.path.abspath(foundFile)

    def _findGeneratedDB(self):
        fname = os.path.join(tmpdir,"gluciferDB"+self._getUniqueName(prefix=self.num)+".gldb")
        if not os.path.isfile(fname):
            raise RuntimeError("The database does not appear to have been created. Please contact developers.")
        
        return os.path.abspath(fname)

    def savefig(self,filename):
        """  Saves the generated image to the provided filename.
            
             Args:
               filename (str):  Filename to save file to.  May include an absolute or relative path.
        """
        self._generateDB()
        self._generateImage()
        generatedFilename=self._findGeneratedFile()
        
        absfilename = os.path.abspath(filename)
        
        # lets set the final extension to that of the glucifer generated file
        splitabsfilename = os.path.splitext(absfilename)
        splitgenfilename = os.path.splitext(generatedFilename)
        if splitabsfilename[1].lower() in [".png", ".jpg", ".jpeg"]:
            frontpart = splitabsfilename[0]
        else:
            frontpart = absfilename
        finaloutFile = frontpart+splitgenfilename[1]
        from shutil import copyfile
        copyfile(generatedFilename,finaloutFile)

    def saveDB(self,filename):
        """  Saves the generated database to the provided filename.
            
             Args:
               filename (str):  Filename to save file to.  May include an absolute or relative path.
        """
        self._generateDB()
        self._generateImage()
        generatedDB=self._findGeneratedDB()
        absfilename = os.path.abspath(filename)

        from shutil import copyfile
        copyfile(generatedDB,absfilename)


    def _generateDB(self):
        vp  = self.componentPointerDictionary[self._localNames[ "vp"]]
        db  = self.componentPointerDictionary[self._localNames[ "db"]]
        win = self.componentPointerDictionary[self._localNames["win"]]
        
        # remove any existing
        for ii in range(vp.drawingObject_Register.objects.count,0,-1):
            _stgermain.StGermain._Stg_ObjectList_RemoveByIndex(vp.drawingObject_Register.objects,ii-1, _stgermain.StGermain.KEEP)
        # first add drawing objects to viewport
        if len(self.drawingObjects) == 0:
            raise RuntimeError("There appears to be no drawing objects to render.")
        for object in self.drawingObjects:
            objectPtr = object._getStgPtr()
            objectPtr.id = 0
            _stgermain.StGermain.Stg_ObjectList_Append(vp.drawingObject_Register.objects,objectPtr)
        # go ahead and fill db
        _stgermain.gLucifer.lucDatabase_DeleteWindows(db)
        _stgermain.gLucifer.lucDatabase_OutputWindow(db, win)
        _stgermain.gLucifer._lucDatabase_Execute(db,None)
        _stgermain.gLucifer._lucWindow_Execute(win,None)

    def _generateImage(self):
        db  = self.componentPointerDictionary[self._localNames[ "db"]]
        # go ahead and draw
        _stgermain.gLucifer.lucDatabase_Dump(db)

    def clear(self):
        """    Clears all the figure's drawing objects 
        """
        del self.drawingObjects[:]

    def __del__(self):
        super(Figure,self).__del__()

    def Surface(self, **kwargs):
        """    Add a surface drawing object to the current figure.
               See 'help(Surface)' for information on the Surface class and it's options.
               
               Returns the generated Surface object.
        """
        guy = _drawing.Surface(**kwargs)
        self.drawingObjects.append(guy)
        return guy

    def Points(self, **kwargs):
        """    Add a points drawing object to the current figure.
               See 'help(Points)' for information on the Points class and it's options.
               
               Returns the generated Points object.
        """
        guy = _drawing.Points(**kwargs)
        self.drawingObjects.append(guy)
        return guy

    def VectorArrows(self, **kwargs):
        """    Add a vector arrow drawing object to the current figure.
               See 'help(VectorArrows)' for information on the VectorArrows class and it's options.
               
               Returns the generated Points object.
        """
        guy = _drawing.VectorArrows(**kwargs)
        self.drawingObjects.append(guy)
        return guy


## weakref dictionary used to allow objects to die
#_figures = weakref.WeakValueDictionary()
_figures = {}

def figure(**kwargs):
    """   Creates instances of Figure object using provided arguments, 
          or returns existing object with provided name.
          
          If no name is provided, an increment integer identifer is generated.
          
          See 'help(Figure)' for information on the Figure class and it's options.
    """
    
    _deleteStale()

    # if name found, check if exists and return
    if "num" in kwargs:
        if kwargs["num"] in _figures:
            return _figures[kwargs["num"]]
    
    # if no figure exists in dict, go ahead and create one
    figGuy = Figure(**kwargs)
    # if id provided, add to dict:
    if "num" in kwargs:
        _figures[figGuy.num] = figGuy

    return figGuy

def close(num=None):
    """ Close figures.  If no argument is provided, close all figures.  
        Note that if there are local references to figures, these figures will remain active.

        Args:
        num (str,int):  identifier of figure to close.  Optional.
    """
    _deleteStale()
    
    if not num:
        for k in _figures.keys():
            del _figures[k]
    else:
        if not num in _figures:
            raise ValueError("Figure with identifier "+num+" not found in list of open figures.")
        del _figures[num]

def _deleteStale():
    # lets first delete any stale items
    for k in _figures.keys():
        if not isinstance(_figures[k],Figure):
            del _figures[k]
