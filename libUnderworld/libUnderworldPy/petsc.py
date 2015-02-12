# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_petsc', [dirname(__file__)])
        except ImportError:
            import _petsc
            return _petsc
        if fp is not None:
            try:
                _mod = imp.load_module('_petsc', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _petsc = swig_import_helper()
    del swig_import_helper
else:
    import _petsc
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        object.__setattr__(self, name, value)
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def OptionsInsertString(string):
    return _petsc.OptionsInsertString(string)
OptionsInsertString = _petsc.OptionsInsertString

def OptionsPrint():
    return _petsc.OptionsPrint()
OptionsPrint = _petsc.OptionsPrint

def OptionsClear():
    return _petsc.OptionsClear()
OptionsClear = _petsc.OptionsClear

def OptionsSetValue(iname, value):
    return _petsc.OptionsSetValue(iname, value)
OptionsSetValue = _petsc.OptionsSetValue

def OptionsClearValue(iname):
    return _petsc.OptionsClearValue(iname)
OptionsClearValue = _petsc.OptionsClearValue

def OptionsInsertFile(file):
    return _petsc.OptionsInsertFile(file)
OptionsInsertFile = _petsc.OptionsInsertFile

def EmptyCall():
    return _petsc.EmptyCall()
EmptyCall = _petsc.EmptyCall
# This file is compatible with both classic and new-style classes.


