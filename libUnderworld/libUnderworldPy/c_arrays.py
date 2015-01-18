# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_c_arrays', [dirname(__file__)])
        except ImportError:
            import _c_arrays
            return _c_arrays
        if fp is not None:
            try:
                _mod = imp.load_module('_c_arrays', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _c_arrays = swig_import_helper()
    del swig_import_helper
else:
    import _c_arrays
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def cdata(*args):
  return _c_arrays.cdata(*args)
cdata = _c_arrays.cdata

def memmove(*args):
  return _c_arrays.memmove(*args)
memmove = _c_arrays.memmove
class DoubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleArray, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _c_arrays.new_DoubleArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _c_arrays.delete_DoubleArray
    __del__ = lambda self : None;
    def __getitem__(self, *args): return _c_arrays.DoubleArray___getitem__(self, *args)
    def __setitem__(self, *args): return _c_arrays.DoubleArray___setitem__(self, *args)
    def cast(self): return _c_arrays.DoubleArray_cast(self)
    __swig_getmethods__["frompointer"] = lambda x: _c_arrays.DoubleArray_frompointer
    if _newclass:frompointer = staticmethod(_c_arrays.DoubleArray_frompointer)
DoubleArray_swigregister = _c_arrays.DoubleArray_swigregister
DoubleArray_swigregister(DoubleArray)

def DoubleArray_frompointer(*args):
  return _c_arrays.DoubleArray_frompointer(*args)
DoubleArray_frompointer = _c_arrays.DoubleArray_frompointer

class FloatArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FloatArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FloatArray, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _c_arrays.new_FloatArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _c_arrays.delete_FloatArray
    __del__ = lambda self : None;
    def __getitem__(self, *args): return _c_arrays.FloatArray___getitem__(self, *args)
    def __setitem__(self, *args): return _c_arrays.FloatArray___setitem__(self, *args)
    def cast(self): return _c_arrays.FloatArray_cast(self)
    __swig_getmethods__["frompointer"] = lambda x: _c_arrays.FloatArray_frompointer
    if _newclass:frompointer = staticmethod(_c_arrays.FloatArray_frompointer)
FloatArray_swigregister = _c_arrays.FloatArray_swigregister
FloatArray_swigregister(FloatArray)

def FloatArray_frompointer(*args):
  return _c_arrays.FloatArray_frompointer(*args)
FloatArray_frompointer = _c_arrays.FloatArray_frompointer

class IntArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntArray, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _c_arrays.new_IntArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _c_arrays.delete_IntArray
    __del__ = lambda self : None;
    def __getitem__(self, *args): return _c_arrays.IntArray___getitem__(self, *args)
    def __setitem__(self, *args): return _c_arrays.IntArray___setitem__(self, *args)
    def cast(self): return _c_arrays.IntArray_cast(self)
    __swig_getmethods__["frompointer"] = lambda x: _c_arrays.IntArray_frompointer
    if _newclass:frompointer = staticmethod(_c_arrays.IntArray_frompointer)
IntArray_swigregister = _c_arrays.IntArray_swigregister
IntArray_swigregister(IntArray)

def IntArray_frompointer(*args):
  return _c_arrays.IntArray_frompointer(*args)
IntArray_frompointer = _c_arrays.IntArray_frompointer

class UnsignedArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UnsignedArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UnsignedArray, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _c_arrays.new_UnsignedArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _c_arrays.delete_UnsignedArray
    __del__ = lambda self : None;
    def __getitem__(self, *args): return _c_arrays.UnsignedArray___getitem__(self, *args)
    def __setitem__(self, *args): return _c_arrays.UnsignedArray___setitem__(self, *args)
    def cast(self): return _c_arrays.UnsignedArray_cast(self)
    __swig_getmethods__["frompointer"] = lambda x: _c_arrays.UnsignedArray_frompointer
    if _newclass:frompointer = staticmethod(_c_arrays.UnsignedArray_frompointer)
UnsignedArray_swigregister = _c_arrays.UnsignedArray_swigregister
UnsignedArray_swigregister(UnsignedArray)

def UnsignedArray_frompointer(*args):
  return _c_arrays.UnsignedArray_frompointer(*args)
UnsignedArray_frompointer = _c_arrays.UnsignedArray_frompointer

# This file is compatible with both classic and new-style classes.


