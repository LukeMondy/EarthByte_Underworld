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
import boundary
import equations
import fields
import geometry
import live
import material
import matrix
import meshing
import physics
import rheology
import shapes
import solvers
import swarms
import utils
import visualisation

# ok, now set this back to the original value
_sys.setdlopenflags( _oldflags )

# add this to allow ctrl+c signal
import signal as _signal
def _signal_handler(_signal, frame):
    print 'You pressed Ctrl+C!'
    sys.exit(0)

_signal.signal(_signal.SIGINT, _signal_handler)


# lets go right ahead and init now.  user can re-init if necessary.
import _stgermain
_stgermain.StgInit()
_rank = _stgermain.getData().rank
_nprocs = _stgermain.getData().nProcs

def rank():
    """
       Returns the rank of the current processors.

       Args:
       None
       Returns:
       rank (unsigned) : Rank of current processor.
       """
    global _rank
    return _rank

def nProcs():
    """
       Returns the number of processors being utilised by the simulation.

       Args:
       None
       Returns:
       nProcs (unsigned) : Number of processors.
       """
    global _nprocs
    return _nprocs
