import underworld
import underworld._stgermain as _stgermain

class Drawing(_stgermain.StgCompoundComponent):
    def __init__(self, num=None, colours="Purple Blue Green Yellow Orange Red".split(), **kwargs):
        if num and not isinstance(num,(str,int,NoneType)):
            raise TypeError("'num' object passed in must be of python type 'str' or 'int'")
        self._num = num

        if not isinstance(colours,(str,list)):
            raise TypeError("'colours' object passed in must be of python type 'str' or 'list'")
        if isinstance(colours,(str)):
            self._colours = colours.split()
        else:
            self._colours = colours

        # build parent
        super(Drawing,self).__init__(**kwargs)

    @property
    def num(self):
        """    num (str,int): integer or string, optional, default: none
            """
        return self._num

    @property
    def colours(self):
        """    colours (list): list of colours to use to draw object.  Should be provided as a list or a string.
        """
        return self._colours

    def _getStgPtr(self):
        """ This function should return the required underlying StGermain object ptr. It may be overwritten by the child class if necessary. """
        return self.componentPointerDictionary[self._localNames["dr"]]

    def _addToStgDict(self):
        # call parents method
        super(Drawing,self)._addToStgDict()

        self._localNames["cm"] = "cm_" + self._getUniqueName(prefix=self.num)
        self._localNames["dr"] = "dr_" + self._getUniqueName(prefix=self.num)

        self.componentDictionary[self._localNames["cm"]] = {
            "Type"          :"lucColourMap",
            "colours"       :" ".join(self.colours),
            "dynamicRange"  :True
        }
        # add an empty(ish) drawing object.  children should fill it out. in particular they need to set 'Type'
        self.componentDictionary[self._localNames["dr"]] = {
            "ColourMap"     :self._localNames["cm"]
        }

    def __del__(self):
        super(Drawing,self).__del__()


class Surface(Drawing):
    """  This drawing object class draws a surface using the provided scalar field.
    """
    def __init__(self, field, **kwargs):
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
        
        # build parent
        super(Surface,self).__init__(**kwargs)

    @property
    def field(self):
        """    field (str): name of live underworld field for which surfaces should be rendered.  Must be a scalar field.
        """
        return self._field

    def _addToStgDict(self):
        # lets build up component dictionary
        # append random string to provided name to ensure unique component names
        # call parents method

        super(Surface,self)._addToStgDict()
        
        drdict = self.componentDictionary[self._localNames["dr"]]
        drdict[         "Type"] = "lucScalarField"
        drdict["FieldVariable"] = self.field
    
    def __del__(self):
        super(Surface,self).__del__()



class Points(Drawing):
    """  This drawing object class draws a swarm of points.
    """
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
                raise ValueError("colourVariable with name '"+colourVariable+"' not found. A live instance must be available before you can create this object.")
            if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StGermain.Variable_Type ) \
               or  _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.SwarmVariable_Type ) :
                raise ValueError("colourVariable with name '"+colourVariable+"' has type '"+ptr.type+"',\n which does not appear to be a child of '"+_stgermain.StGermain.Variable_Type+"' or '"+_stgermain.StgDomain.SwarmVariable_Type+"' types, as required.")

        if sizeVariable:
            if not isinstance(sizeVariable,(str)):
                raise TypeError("'sizeVariable' object passed in must be of python type 'str'")
            ptr = _stgermain.GetLiveComponent(sizeVariable)
            if not ptr:
                raise ValueError("sizeVariable with name '"+sizeVariable+"' not found. A live instance must be available before you can create this object.")
            if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StGermain.Variable_Type ) \
               or  _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.SwarmVariable_Type ) :
                raise ValueError("sizeVariable with name '"+sizeVariable+"' has type '"+ptr.type+"',\n which does not appear to be a child of '"+_stgermain.StGermain.Variable_Type+"' or '"+_stgermain.StgDomain.SwarmVariable_Type+"' types, as required.")

        if opacityVariable:
            if not isinstance(opacityVariable,(str)):
                raise TypeError("'opacityVariable' object passed in must be of python type 'str'")
            ptr = _stgermain.GetLiveComponent(opacityVariable)
            if not ptr:
                raise ValueError("opacityVariable with name '"+opacityVariable+"' not found. A live instance must be available before you can create this object.")
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
    
    def _addToStgDict(self):
        # lets build up component dictionary

        # call parents method
        super(Points,self)._addToStgDict()

        drdict = self.componentDictionary[self._localNames["dr"]]
        drdict[           "Type"] = "lucSwarmViewer"
        drdict[          "Swarm"] = self.swarm
        drdict[ "ColourVariable"] = self.colourVariable
        drdict[   "SizeVariable"] = self.sizeVariable
        drdict["opacityVariable"] = self.opacityVariable
        drdict[      "pointSize"] = self.pointSize
    

    def __del__(self):
        super(Points,self).__del__()


class VectorArrows(Drawing):
    """  This drawing object class draws vector arrows corresponding to the provided vector field.
    """
    def __init__(self, vectorField, resolutionX=16, resolutionY=16, resolutionZ=16, arrowHeadSize=0.3, lengthScale=0.3,glyphs=3, **kwargs):
        if not isinstance(vectorField,(str)):
            raise TypeError("'vectorField' object passed in must be of python type 'str'")
        ptr = _stgermain.GetLiveComponent(vectorField)
        if not ptr:
            raise ValueError("FieldVariable with name '"+vectorField+"' not found. A live instance must be available before you can create this object.")
        if not _stgermain.StGermain.Stg_Class_CompareType( ptr, _stgermain.StgDomain.FieldVariable_Type ):
            raise ValueError("FieldVariable with name '"+vectorField+"' not a child of '"+_stgermain.StgDomain.FieldVariable_Type+"' type.")
        if not ptr.fieldComponentCount > 1:
            raise ValueError("FieldVariable with name '"+vectorField+"' is not a vector field. It appears to have "+str(ptr.fieldComponentCount)+" components.")

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

        self._vectorField = vectorField
        self._resolutionX = resolutionX
        self._resolutionY = resolutionY
        self._resolutionZ = resolutionZ
        self._arrowHeadSize = arrowHeadSize
        self._lengthScale = lengthScale
        self._glyphs = glyphs

        # build parent
        super(VectorArrows,self).__init__(**kwargs)
    
    @property
    def vectorField(self):
        """    vectorField (str): name of live underworld vector field for which vector arrows will be rendered.
        """
        return self._vectorField
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
    
    def _addToStgDict(self):
        # lets build up component dictionary

        # call parents method
        super(VectorArrows,self)._addToStgDict()

        drdict = self.componentDictionary[self._localNames["dr"]]
        drdict[           "Type"] = "lucVectorArrows"
        drdict[  "FieldVariable"] = self.vectorField
        drdict[    "resolutionX"] = self.resolutionX
        drdict[    "resolutionY"] = self.resolutionY
        drdict[    "resolutionZ"] = self.resolutionZ
        drdict[  "arrowHeadSize"] = self.arrowHeadSize
        drdict[    "lengthScale"] = self.lengthScale
        drdict[         "glyphs"] = self.glyphs
    

    def __del__(self):
        super(VectorArrows,self).__del__()
