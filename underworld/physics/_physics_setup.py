# Physics package
import underworld as _uw

##############################################################################
## This code adds what is required to the python dictionary 
## to configure physics
## We eventually pass the python dictionary back to Underworld
## and Underworld then uses this information to configure and set
## itself up.
##############################################################################

'''
This code adds what is required to the python dictionary
Ultimately the global Dictionary gets passed back to Underworld which then actually creates the simulation

'''

# This is 'physics' but needs a force Vector and an integration Swarm
# These are things we would like to hide from casual users.
# So this function would have to be smarter and find the vector and swarm without users
# having to supply them



def addBuoyancy(forceVector="mom_force", intSwarm="", temperatureField="TemperatureField"):

    globalDict = _uw.GetCurrentPythonDictionary()

    
    buoyancy = _uw.NewComponentEntryInStgDict( globalDict, 
                                               name = "buoyancyForceTerm",
                                               Type = "BuoyancyForceTerm",
                                               TemperatureField = temperatureField,  #optional: temp of 0.0 used if no Field
                                               Swarm = intSwarm,
                                               gravity = "gravity"
                                               )
    return buoyancy
