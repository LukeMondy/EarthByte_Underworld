# borrowed from http://stackoverflow.com/questions/1389180/python-automatically-initialize-instance-variables

from functools import wraps
import inspect

def initializer(func):
    """
    Automatically assigns the parameters.

    >>> class process:
    ...     @initializer
    ...     def __init__(self, cmd, reachable=False, user='root'):
    ...         pass
    >>> p = process('halt', True)
    >>> p.cmd, p.reachable, p.user
    ('halt', True, 'root')
    """
    names, varargs, keywords, defaults = inspect.getargspec(func)

    @wraps(func)
    def wrapper(self, *args, **kargs):
        for name, arg in list(zip(names[1:], args)) + list(kargs.items()):
            setattr(self, name, arg)

        for name, default in zip(reversed(names), reversed(defaults)):
            if not hasattr(self, name):
                setattr(self, name, default)

        func(self, *args, **kargs)

    return wrapper

def pointerToArray( swigObject, shape, type=None):
    """"
    Casts a swig pointer into a numpy array.  Only native types (int, float, etc) are supported.
    Note that no copying occurs here, and all values may be modified in place.
    
    Arguments:
        swigObject (swig type): pointer to the underlying object you wish to cast as a numpy array.
                                         note that type must be native, such as int or float.
        shape (tuple):  Shape of underlying array, provided as a tuple.
                             Note that providing an incorrect shape will allow users to modify memory out of bounds
                             of the underlying c allocation, which is very naughty.
        type (str), optional:  Explicit type to use.  Options are 'char', 'uint', 'int', 'float', 'double'.
                             
    Returns:
        numpy array (ndarray):  Required numpy array.
        
    """
    import numpy as np
    import ctypes

    if type == None:
        # lets try and find it's type
        typestr = repr(swigObject).split("'")[1]

        if typestr == 'char *':
            type = 'char'
        elif typestr == 'unsigned int *':
            type = 'uint'
        elif typestr == 'int *':
            type = 'int'
        elif typestr == 'float *':
            type = 'float'
        elif typestr == 'double *':
            type = 'double'
        else:
            raise ValueError("The swig pointer you passed in appears to have type '{}' which is not supported.\n"
                              "If you believe this is a supported type, please explicitly provide a type.".format(typestr) )

    if type not in ['char', 'uint', 'int', 'float', 'double']:
        raise ValueError( "Provided type '{}' not recognised.\nAvailable types are 'char', 'uint', 'int', 'float' and 'double'".format(type))

    try:
        longPtr = swigObject.this.__long__()
    except:
        try:
            longPtr = swigObject.__long__()
        except:
            raise ValueError("Unable to extract pointer from provided swig object.  If you believe this is an error, please contact developers.")



    if type == 'char':
        cptr = ctypes.cast( longPtr, ctypes.POINTER( ctypes.c_char ) )
    elif type == 'uint':
        cptr = ctypes.cast( longPtr, ctypes.POINTER( ctypes.c_uint ) )
    elif type == 'int':
        cptr = ctypes.cast( longPtr, ctypes.POINTER( ctypes.c_int ) )
    elif type == 'float':
        cptr = ctypes.cast( longPtr, ctypes.POINTER( ctypes.c_float ) )
    elif type == 'double':
        cptr = ctypes.cast( longPtr, ctypes.POINTER( ctypes.c_double ) )
    else:
        raise RuntimeError("Something has gone wrong. Please contact developers." )

    return np.ctypeslib.as_array(cptr, shape=shape)

