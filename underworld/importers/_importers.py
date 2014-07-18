import underworld
import underworld._stgermain as _stgermain
import numpy as np

class VoxelDataHandler_ndarray(_stgermain.StgCompoundComponent):
    """
    This Class wraps the VoxelDataHandler_ndarray StGermain class.
    This can be used to view multidimensional numpy arrays as voxel datasets within Underworld.
    """
    def __init__(self, ndarray, minTup=(0.,0.,0.), maxTup=(1.,1.,1.), **kwargs):
        if not isinstance(ndarray,(np.ndarray)):
            raise TypeError("'ndarray' object passed in must be of type 'ndarray'")
        if len(ndarray.shape) < 2 or len(ndarray.shape) > 3:
            raise ValueError("'ndarray' object must be of dimensionality 2 or 3.")
        self._ndarray = ndarray

        if not isinstance(minTup,(tuple)):
            raise TypeError("'minTup' object passed in must be of python type 'tuple'")
        if len(minTup) < 2 or len(minTup) > 3:
            raise ValueError("'minTup' object must be of dimensionality 2 or 3.")
        for datum in minTup:
            if not isinstance(datum,(int,float)):
                raise TypeError("'minTup' object passed in must only contain python 'int' or 'float' data")
        self._minTup = minTup

        if not isinstance(maxTup,(tuple)):
            raise TypeError("'maxTup' object passed in must be of python type 'tuple'")
        if len(maxTup) < 2 or len(maxTup) > 3:
            raise ValueError("'maxTup' object must be of dimensionality 2 or 3.")
        for datum in maxTup:
            if not isinstance(datum,(int,float)):
                raise TypeError("'maxTup' object passed in must only contain python 'int' or 'float' data")
        self._maxTup = maxTup

        # build parent
        super(VoxelDataHandler_ndarray,self).__init__(**kwargs)

    @property
    def ndarray(self):
        """    ndarray (ndarray): numpy array for VoxelDataHandler_ndarray to utilise
        """
        return self._ndarray

    @property
    def minTup(self):
        """    minTup (tuple(int,float)): location to be considered the lower bound on the dataset domain.
        """
        return self._minTup
    @property
    def maxTup(self):
        """    maxTup (tuple(int,float)): location to be considered the upper bound on the dataset domain.
        """
        return self._maxTup

    def _addToStgDict(self):
        # call parents method
        super(VoxelDataHandler_ndarray,self)._addToStgDict()

        self._localNames["vdh_np"] = "vdh_np_" + self._getUniqueName()

        if len(self.ndarray.shape) == 2:
            numcellsk = 1
        else:
            numcellsk = self.ndarray.shape[2]

        if len(self.minTup) == 2:
            startK = 0
        else:
            startK = self.minTup[2]

        if len(self.maxTup) == 2:
            finK = 1
        else:
            finK = self.maxTup[2]

        cellSizeI = float(self.maxTup[0] - self.minTup[0])/float(self.ndarray.shape[0])
        cellSizeJ = float(self.maxTup[1] - self.minTup[1])/float(self.ndarray.shape[1])
        cellSizeK = float(finK - startK)/float(numcellsk)

        if   np.issubdtype(self.ndarray.dtype,np.int8):
            nptype = "char"
        elif np.issubdtype(self.ndarray.dtype,np.int32):
            nptype = "int"
        elif np.issubdtype(self.ndarray.dtype,np.float32):
            nptype = "float"
        elif np.issubdtype(self.ndarray.dtype,np.float64):
            nptype = "double"
        else:
            raise ValueError("Provided numpy array does not appear to be of a supported type.\n"+\
                             "Type is "+str(self.ndarray.dtype)+" while supported types are 'int8', 'int32', 'float32' and 'float64'")

        self.componentDictionary[self._localNames["vdh_np"]] = {
            "Type"              :"VoxelDataHandler_ndarray",
            "ndPointer"         :hex(self.ndarray.__array_interface__['data'][0]),  # note we convert the pointer to ndarray data to a string here
            "NumCellsI"         :self.ndarray.shape[0],
            "NumCellsJ"         :self.ndarray.shape[1],
            "NumCellsK"         :numcellsk,
            "StartCoordI"       :self.minTup[0]+0.5*cellSizeI,
            "StartCoordJ"       :self.minTup[1]+0.5*cellSizeJ,
            "StartCoordK"       :startK+0.5*cellSizeK,
            "CellSizeI"         :cellSizeI,
            "CellSizeJ"         :cellSizeJ,
            "CellSizeK"         :cellSizeK,
            "DataType"          :nptype,
            "mapIAxisToStgAxis" :"X",
            "mapJAxisToStgAxis" :"Y",
            "mapKAxisToStgAxis" :"Z"
        }
        super(VoxelDataHandler_ndarray,self)._addToStgDict()

    def __del__(self):
        super(VoxelDataHandler_ndarray,self).__del__()


class NumpyVoxelField(_stgermain.StgCompoundComponent):
    """
    This Class takes a numpy multidimensional array, and makes it available as an Underworld Field.
    """
    __doc__ += VoxelDataHandler_ndarray.__doc__  # this adds docstring info from VoxelDataHandler_ndarray class
    def __init__(self, *args, **kwargs):
    
        self.vdh_np = VoxelDataHandler_ndarray(*args, **kwargs)
        
        # build parent
        super(NumpyVoxelField,self).__init__(**kwargs)


    def _addToStgDict(self):
        # call parents method
        super(NumpyVoxelField,self)._addToStgDict()

        self._localNames["nvf"] = "nvf" + self._getUniqueName()

        self.componentDictionary[self._localNames["nvf"]] = {
            "Type"                      :"VoxelFieldVariable",
            "VoxelDataHandler"          :self.vdh_np._localNames["vdh_np"],
            "UseNearestCellIfOutside"   :True,
            "fieldComponentCount"       :1,
            "dim"                       :len(self.vdh_np.ndarray.shape) # get the numpy shape for the dimensionality
        }

    def __del__(self):
        super(NumpyVoxelField,self).__del__()
