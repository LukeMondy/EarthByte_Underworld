
# Solvers module - activates and configures the "new solver system"

'''
    This would be activated by the following XML / opt files

    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/VelocityMassMatrixSLE.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
    -options_file $UWPATH/Solvers/examples/options-scr-mg-accelerating.opt \

    We try to provide a more flexible interface to this lot

'''
import setup
import live 
## These are all just dictionary setup as well 
# from _multigrid import *
#from _options import *

# Initialisation

_solversOn = False  # obsolete
_solverType = "None"

def solverType():
    """
        Returns a string detailing the solver
    """
    return _solverType

def setSolverType( typeString ):
    """
    Sets the string representing the solver

    """    
    
    global _solverType
    _solverType = typeString

    return    
