import underworld
from underworld import initializer
import underworld._stgermain as _stgermain
import numpy as np

class Drawing(_stgermain.StgCompoundComponent):
    def __new__(cls, objectDict, *args, **kwargs):
        
        if not isinstance(objectDict, dict):
            raise TypeError("objectDict passed in must be of python type 'dict' or subclass")
        
        if "dr" not in objectDict:
            raise ValueError("objectDict passed in from child must contain an object with key 'dr' mapping to a value which defines the StGermain object type.")
        
        # also add colourmap object
        objectDict["cm"] = "lucColourMap"

        return super(Drawing,cls).__new__(cls, objectDict, *args, **kwargs)
    
    def __init__(self, num=None, colours="Purple Blue Green Yellow Orange Red".split(), opacity=-1, **kwargs):
        self.num = num
        self.colours = colours
        self.opacity = opacity
        
        # build parent
        super(Drawing,self).__init__()

    def _addToStgDict(self):
        
        self.componentDictionary[self._clib_cm.name] = {
            "colours"       :" ".join(self.colours),
            "dynamicRange"  :True
        }
        # add an empty(ish) drawing object.  children should fill it out.
        self.componentDictionary[self._clib_dr.name] = {
            "ColourMap"     :self._clib_cm.name,
            "opacity"       :self.opacity
        }

        # call parents method
        super(Drawing,self)._addToStgDict()

    def __del__(self):
        super(Drawing,self).__del__()


    @property
    def num(self):
        """    num (str,int): integer or string identifier. Must not contain spaces. optional, default: none
        """
        return self._num
    @num.setter
    def num(self, num):
        self._setterAssertNotConstructed()
        if num and not isinstance(num,(str,int,NoneType)):
            raise TypeError("'num' object passed in must be of python type 'str' or 'int'")
        if num and isinstance(num,str) and (" " in num):
            raise ValueError("'num' object passed in must not contain any spaces.")
        self._num = num

    @property
    def colours(self):
        """    colours (list): list of colours to use to draw object.  Should be provided as a list or a string.
        """
        return self._colours
    @colours.setter
    def colours(self,colours):
        self._setterAssertNotConstructed()
        if not isinstance(colours,(str,list)):
            raise TypeError("'colours' object passed in must be of python type 'str' or 'list'")
        if isinstance(colours,(str)):
            self._colours = colours.split()
        else:
            self._colours = colours

    @property
    def opacity(self):
        """    opacity (float): Opacity of drawing object.  Takes values from 0. to 1., while a value of -1 explicitly disables opacity.
        """
        return self._opacity
    @opacity.setter
    def opacity(self,opacity):
        self._setterAssertNotConstructed()
        if not isinstance(opacity,(int,float)):
            raise TypeError("'opacity' object passed in must be of python type 'int' or 'float'")
        if float(opacity) > 1. or float(opacity) < -1.:
            raise ValueError("'opacity' object must takes values from 0. to 1., while a value of -1 explicitly disables opacity.")
        self._opacity = opacity

class Surface(Drawing):
    """  This drawing object class draws a surface using the provided scalar field.
    """
    def __new__(cls, objectDict={}, *args, **kwargs):
        
        if not isinstance(objectDict, dict):
            raise TypeError("objectDict passed in must be of python type 'dict' or subclass")
        
        # note we check if it already exists incase a child object has set it
        if "dr" not in objectDict:
            objectDict["dr"] = "lucScalarField"
        
        return super(Surface,cls).__new__(cls, objectDict, *args, **kwargs)

    def __init__(self, field=None, ndarray=None, drawSides="xyzXYZ",*args, **kwargs):
        if not((field==None) or (ndarray==None)):
            raise ValueError("Either an Underworld field or a numpy array must be set as arguments when initialising a Surface object, but not both.")

        if field:
            if not isinstance(field,(str)):
                raise TypeError("'field' object passed in must be of python type 'str'")
            fieldPtr = _stgermain.GetLiveComponent(field)
            if not fieldPtr:
                raise ValueError("Field with name '"+field+"' not found. A live instance must be available before you can create this object.")
            if not _stgermain.StGermain.Stg_Class_CompareType( fieldPtr, _stgermain.StgDomain.FieldVariable_Type ):
                raise ValueError("Field with name '"+field+"' not a child of '"+_stgermain.StgDomain.FieldVariable_Type+"' type.")
            if not fieldPtr.fieldComponentCount == 1:
                raise ValueError("Field with name '"+field+"' is not a scalar field. It appears to have "+str(fieldPtr.fieldComponentCount)+" components.")
        self._field = field
        
        if not ndarray==None:
            self._nvf = underworld.importers.NumpyVoxelField(ndarray=ndarray, **kwargs)
        self._ndarray=ndarray
        
        if not isinstance(drawSides,str):
            raise ValueError("'drawSides' argument must be of python type 'str'")
        self._drawSides = drawSides
        
        # build parent
        super(Surface,self).__init__(**kwargs)

    def _addToStgDict(self):
        # lets build up component dictionary
        # append random string to provided name to ensure unique component names
        # call parents method
        
        drdict = self.componentDictionary[self._clib_dr.name]
        drdict["drawSides"] = self.drawSides
        if self.field:
            drdict["FieldVariable"] = self.field
        else:
            drdict["FieldVariable"] = self._nvf._clib_nvf.name

        super(Surface,self)._addToStgDict()

    def __del__(self):
        super(Surface,self).__del__()


    @property
    def field(self):
        """    field (str): name of live underworld field for which surfaces should be rendered.  Must be a scalar field.
        """
        return self._field

    @property
    def ndarray(self):
        """    ndarray (numpy.ndarray): numpy ndarray object for which surfaces should be rendered.  Must be of dimensionality 2 or 3.
        """
        return self._ndarray

    @property
    def drawSides(self):
        """    drawSides (str): sides (x,y,z,X,Y,Z) for which the surface should be drawn.  default is all sides ("xyzXYZ").
        """
        return self._drawSides


class Points(Drawing):
    """  This drawing object class draws a swarm of points.
    """
    def __new__(cls, objectDict={}, *args, **kwargs):
        
        if not isinstance(objectDict, dict):
            raise TypeError("objectDict passed in must be of python type 'dict' or subclass")
        
        # note we check if it already exists incase a child object has set it
        if "dr" not in objectDict:
            objectDict["dr"] = "lucSwarmViewer"
        
        return super(Points,cls).__new__(cls, objectDict, *args, **kwargs)

    def __init__(self, swarm, colourVariable=None, sizeVariable=None, opacityVariable=None, pointSize=1.0, **kwargs):
        if not isinstance(swarm,(str)):
            raise TypeError("'swarm' object passed in must be of python type 'str'")
        ptr = _stgermain.GetLiveComponent(swarm)
        if not ptr:
            raise ValueError("Swarm with name '"+swarm+"' not found. A live instance must be available before you can create this object.")
        if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.Swarm_Type ):
            raise ValueError("Swarm with name '"+swarm+"' not a child of '"+_stgermain.StgDomain.Swarm_Type+"' type.")

        if colourVariable:
            if not isinstance(colourVariable,(str)):
                raise TypeError("'colourVariable' object passed in must be of python type 'str'")
            ptr = _stgermain.GetLiveComponent(colourVariable)
            if not ptr:
                raise ValueError("colourVariable with name '"+colourVariable+"' not found. A live instance must be available before you can create this object.\n"\
                                 "Run 'underworld.swarms.tools.Swarm_PrintVariables(\""+swarm+"\")' to find available swarm variables.")
            if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StGermain.Variable_Type ) \
               or  _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.SwarmVariable_Type ) :
                raise ValueError("colourVariable with name '"+colourVariable+"' has type '"+ptr.type+"',\n which does not appear to be a child of '"+_stgermain.StGermain.Variable_Type+"' or '"+_stgermain.StgDomain.SwarmVariable_Type+"' types, as required.")

        if sizeVariable:
            if not isinstance(sizeVariable,(str)):
                raise TypeError("'sizeVariable' object passed in must be of python type 'str'")
            ptr = _stgermain.GetLiveComponent(sizeVariable)
            if not ptr:
                raise ValueError("sizeVariable with name '"+sizeVariable+"' not found. A live instance must be available before you can create this object.\n"\
                                 "Run 'underworld.swarms.tools.Swarm_PrintVariables(\""+swarm+"\")' to find available swarm variables.")
            if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StGermain.Variable_Type ) \
               or  _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.SwarmVariable_Type ) :
                raise ValueError("sizeVariable with name '"+sizeVariable+"' has type '"+ptr.type+"',\n which does not appear to be a child of '"+_stgermain.StGermain.Variable_Type+"' or '"+_stgermain.StgDomain.SwarmVariable_Type+"' types, as required.")

        if opacityVariable:
            if not isinstance(opacityVariable,(str)):
                raise TypeError("'opacityVariable' object passed in must be of python type 'str'")
            ptr = _stgermain.GetLiveComponent(opacityVariable)
            if not ptr:
                raise ValueError("opacityVariable with name '"+opacityVariable+"' not found. A live instance must be available before you can create this object.\n"\
                                 "Run 'underworld.swarms.tools.Swarm_PrintVariables(\""+swarm+"\")' to find available swarm variables.")
            if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StGermain.Variable_Type ) \
               or  _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.SwarmVariable_Type ) :
                raise ValueError("opacityVariable with name '"+opacityVariable+"' has type '"+ptr.type+"',\n which does not appear to be a child of '"+_stgermain.StGermain.Variable_Type+"' or '"+_stgermain.StgDomain.SwarmVariable_Type+"' types, as required.")

        if not isinstance(pointSize,(float,int)):
            raise TypeError("'pointSize' object passed in must be of python type 'float'")

        self._swarm = swarm
        self._colourVariable = colourVariable
        self._sizeVariable = sizeVariable
        self._opacityVariable = opacityVariable
        self._pointSize = pointSize

        # build parent
        super(Points,self).__init__(**kwargs)

    def _addToStgDict(self):
        # lets build up component dictionary
        
        # call parents method
        super(Points,self)._addToStgDict()
        
        drdict = self.componentDictionary[self._clib_dr.name]
        drdict[          "Swarm"] = self.swarm
        drdict[ "ColourVariable"] = self.colourVariable
        drdict[   "SizeVariable"] = self.sizeVariable
        drdict["opacityVariable"] = self.opacityVariable
        drdict[      "pointSize"] = self.pointSize


    def __del__(self):
        super(Points,self).__del__()

    @property
    def swarm(self):
        """    swarm (str): name of live underworld swarm for which points will be rendered.
        """
        return self._swarm
    @property
    def colourVariable(self):
        """    colourVariable (str): name of live underworld swarm variable which will determine the point colours.
        """
        return self._colourVariable
    @property
    def sizeVariable(self):
        """    sizeVariable (str): name of live underworld swarm variable which will determine the point size.
        """
        return self._sizeVariable
    @property
    def opacityVariable(self):
        """    opacityVariable (str): name of live underworld swarm variable which will determine the point opacity.
        """
        return self._opacityVariable

    @property
    def pointSize(self):
        """    pointSize (float): size of points
        """
        return self._pointSize


class VectorArrows(Drawing):
    """  This drawing object class draws vector arrows corresponding to the provided vector field.
    """
    
    def __new__(cls, objectDict={}, *args, **kwargs):
        
        if not isinstance(objectDict, dict):
            raise TypeError("objectDict passed in must be of python type 'dict' or subclass")
        
        # note we check if it already exists incase a child object has set it
        if "dr" not in objectDict:
            objectDict["dr"] = "lucVectorArrows"
        
        return super(VectorArrows,cls).__new__(cls, objectDict, *args, **kwargs)


    def __init__(self, field, resolutionX=16, resolutionY=16, resolutionZ=16, arrowHeadSize=0.3, lengthScale=0.3,glyphs=3, **kwargs):
        if not isinstance(field,(str)):
            raise TypeError("'field' object passed in must be of python type 'str'")
        ptr = _stgermain.GetLiveComponent(field)
        if not ptr:
            raise ValueError("FieldVariable with name '"+field+"' not found. A live instance must be available before you can create this object.")
        if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.FieldVariable_Type ):
            raise ValueError("FieldVariable with name '"+field+"' not a child of '"+_stgermain.StgDomain.FieldVariable_Type+"' type.")
        if not ptr.fieldComponentCount > 1:
            raise ValueError("FieldVariable with name '"+field+"' is not a vector field. It appears to have "+str(ptr.fieldComponentCount)+" components.")

        if resolutionX:
            if not isinstance(resolutionX,(int)):
                raise TypeError("'resolutionX' object passed in must be of python type 'int'")
        if resolutionY:
            if not isinstance(resolutionY,(int)):
                raise TypeError("'resolutionY' object passed in must be of python type 'int'")
        if resolutionZ:
            if not isinstance(resolutionZ,(int)):
                raise TypeError("'resolutionZ' object passed in must be of python type 'int'")
        if arrowHeadSize:
            if not isinstance(arrowHeadSize,(float,int)):
                raise TypeError("'arrowHeadSize' object passed in must be of python type 'int' or 'float'")
            if arrowHeadSize < 0 or arrowHeadSize > 1:
                raise ValueError("'arrowHeadSize' can only take values between zero and one. Value provided is " + str(arrowHeadSize)+".")
        if lengthScale:
            if not isinstance(lengthScale,(float,int)):
                raise TypeError("'lengthScale' object passed in must be of python type 'int' or 'float'")
        if glyphs:
            if not isinstance(glyphs,(int)):
                raise TypeError("'glyphs' object passed in must be of python type 'int'")

        self._field = field
        self._resolutionX = resolutionX
        self._resolutionY = resolutionY
        self._resolutionZ = resolutionZ
        self._arrowHeadSize = arrowHeadSize
        self._lengthScale = lengthScale
        self._glyphs = glyphs

        # build parent
        super(VectorArrows,self).__init__(**kwargs)

    def _addToStgDict(self):
        # lets build up component dictionary
        
        # call parents method
        super(VectorArrows,self)._addToStgDict()
        
        drdict = self.componentDictionary[self._clib_dr.name]
        drdict[  "FieldVariable"] = self.field
        drdict[    "resolutionX"] = self.resolutionX
        drdict[    "resolutionY"] = self.resolutionY
        drdict[    "resolutionZ"] = self.resolutionZ
        drdict[  "arrowHeadSize"] = self.arrowHeadSize
        drdict[    "lengthScale"] = self.lengthScale
        drdict[         "glyphs"] = self.glyphs


    def __del__(self):
        super(VectorArrows,self).__del__()


    @property
    def field(self):
        """    field (str): name of live underworld vector field for which vector arrows will be rendered.
        """
        return self._field
    @property
    def resolutionX(self):
        """    resolutionX (int): Number of vector arrows to render in the X direction. Default is 16.
        """
        return self._resolutionX
    @property
    def resolutionY(self):
        """    resolutionY (int): Number of vector arrows to render in the Y direction. Default is 16.
        """
        return self._resolutionY
    @property
    def resolutionZ(self):
        """    resolutionZ (int): Number of vector arrows to render in the Z direction. Default is 16.
        """
        return self._resolutionZ
    @property
    def arrowHeadSize(self):
        """    arrowHeadSize (float): The size of the head of the arrow compared with the arrow length. Must be between [0, 1].   Default is 0.3.
        """
        return self._arrowHeadSize
    @property
    def lengthScale(self):
        """    lengthScale (float): A factor to scale the size of the arrows by.  Default is 0.3.
        """
        return self._lengthScale
    @property
    def glyphs(self):
        """    glyphs (int): Type of glyph to render for vector arrow.
        """
        return self._glyphs
