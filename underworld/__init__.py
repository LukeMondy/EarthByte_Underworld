# ok, first need to change default python dlopen flags to global
# this is because when python imports the module, the shared libraries are loaded as RTLD_LOCAL
# and then when MPI_Init is called, OpenMPI tries to dlopen its plugin, they are unable to 
# link to the openmpi libraries as they are private
import sys as _sys
import ctypes as _ctypes
_oldflags = _sys.getdlopenflags()
_sys.setdlopenflags( _oldflags | _ctypes.RTLD_GLOBAL )

import libUnderworld

from _underworld import *

import visualisation

# ok, now set this back to the original value
_sys.setdlopenflags( _oldflags )

# add this to allow ctrl+c signal
import signal as _signal
def _signal_handler(_signal, frame):
        print 'You pressed Ctrl+C!'
        sys.exit(0)

_signal.signal(_signal.SIGINT, _signal_handler)



