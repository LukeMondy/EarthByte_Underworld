#!/usr/bin/env python
'''
In this example it is show how you might perform a parametric simulation,
keeping Underworld alive between individual solves.  This avoids the requirement
of tearing down and rebuilding the simulation (particles/meshes/etc) between every solve,
which potentially will save significant time.

Note in this example we directly interact with c-struct objects to change in this case the boundary value.
This is a more advanced usage pattern as it is quiet easy to make erroneous modifications.  Where a complete teardown
does not impose significant costs, it is recommended that users instead effect changes through the input dictionary.

'''

import uwpytools


# init using json
uwpytools.InitWithArgs("SimpleThermal.xml")

stgdict = uwpytools.GetCurrentDictionary()

uwpytools.Construct()
uwpytools.BuildAndInitialise()

# grab the context
#context = uwpytools.GetLiveComponent("context")

# grab the LowerBoundaryPpc
LowerBoundaryPpc = uwpytools.GetLiveComponent("LowerBoundaryPpc")

# grab the TemperatureField
TemperatureField = uwpytools.GetLiveComponent("TemperatureField")

temperatureMesh = TemperatureField.feMesh


print "Temperature field - ", TemperatureField
print "Temperature mesh  - ", temperatureMesh
print "MinSep            - ", temperatureMesh.minSep
print "Temperature verts - ", temperatureMesh.verts  
print "Temperature vars  - ", temperatureMesh.vars    



# lets create some values to set as the BC
lowerBCVals = [3000,300,30,3,.3]

for val in lowerBCVals:
	# note that here we directly modify the Ppc c struct value 'value'.
	LowerBoundaryPpc.value = val
	# do the solve
	uwpytools.Step()
	# print some things
	if uwpytools.rank() == 0:
		print "------------------------------------"
		print "TemperatureField min value", uwpytools.FieldVariable_GetMinFieldMagnitude(TemperatureField)
		print "TemperatureField max value", uwpytools.FieldVariable_GetMaxFieldMagnitude(TemperatureField)
		print "------------------------------------"


uwpytools.Finalise()
