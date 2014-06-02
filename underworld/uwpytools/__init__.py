# ok, first need to change default python dlopen flags to global
# this is because when python imports the module, the shared libraries are loaded as RTLD_LOCAL
# and then when MPI_Init is called, OpenMPI tries to dlopen its plugin, they are unable to 
# link to the openmpi libraries as they are private
import sys as _sys
import ctypes as _ctypes
_oldflags = _sys.getdlopenflags()
_sys.setdlopenflags( _oldflags | _ctypes.RTLD_GLOBAL )


### This allows us to see the public contents of _uwpytools 
### directly under the uwpytools namespace. If you need access
### to the private variables of that file for any submodule, 
### you can access them through uwpytools._uwpytools.NAME 

from _uwpytools import *

import StGermain
import StgDomain
import StgFEM
import PICellerator
import Underworld
import gLucifer
import c_arrays
import c_pointers

import visualisation

# ok, now set this back to the original value
_sys.setdlopenflags( _oldflags )

# add this to allow ctrl+c signal
import signal as _signal
def _signal_handler(_signal, frame):
        print 'You pressed Ctrl+C!'
        sys.exit(0)

_signal.signal(_signal.SIGINT, _signal_handler)



