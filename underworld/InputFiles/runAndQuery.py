#!/usr/bin/env python

'''
This example shows how you might control the simulation flow, obtaining various metrics at certain timesteps.  
You might also modify various aspects of the simulation at different times in the simulation. 

'''



import uwpytools

# init using json
uwpytools.InitWithArgs("RayleighTaylorBenchmark.json")

uwpytools.Construct()
uwpytools.BuildAndInitialise()

# grab the context
context = uwpytools.GetLiveComponent("context")
# grab the velMag field
VelocityMagnitudeField = uwpytools.GetLiveComponent("VelocityMagnitudeField")

minRange,maxRange = uwpytools.FieldVariable_GetMinAndMaxLocalCoords(VelocityMagnitudeField)

while context.timeStep < 2:
	uwpytools.Step()


print ""
print ""
print ""
if uwpytools.rank() == 0:
	print "------------------------------------"
	print "VelocityMagnitudeField local domain (min,max)", uwpytools.FieldVariable_GetMinAndMaxLocalCoords(VelocityMagnitudeField)
	print "VelocityMagnitudeField global domain (min,max)", uwpytools.FieldVariable_GetMinAndMaxGlobalCoords(VelocityMagnitudeField)
	print "------------------------------------"
print ""
print ""
# do some things
if uwpytools.rank() == 0:
	print "------------------------------------"
	print "At timestep", context.timeStep
	print "VelocityMagnitudeField min value", uwpytools.FieldVariable_GetMinFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField max value", uwpytools.FieldVariable_GetMaxFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField Integral", uwpytools.FeVariable_Integrate(VelocityMagnitudeField)

	print "VelocityMagnitudeField interpolation results:"
	for k in range(0,10):
		xcoord = minRange[0] + float(k)/10 * (maxRange[0]-minRange[0])
		ycoord = minRange[1] + float(k)/10 * (maxRange[1]-minRange[1])
		print "    coord=", (xcoord,ycoord), "result=", uwpytools.FieldVariable_InterpolateValueAt(VelocityMagnitudeField, (xcoord,ycoord) )
	print "------------------------------------"
print ""
print ""
print ""

# run for a few more steps
while context.timeStep < 4:
	uwpytools.Step()

print ""
print ""
print ""
if uwpytools.rank() == 0:
	print "------------------------------------"
	print "At timestep", context.timeStep
	print "VelocityMagnitudeField min value", uwpytools.FieldVariable_GetMinFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField max value", uwpytools.FieldVariable_GetMaxFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField Integral", uwpytools.FeVariable_Integrate(VelocityMagnitudeField)

	print "VelocityMagnitudeField interpolation results:"
	for k in range(0,10):
		xcoord = minRange[0] + float(k)/10 * (maxRange[0]-minRange[0])
		ycoord = minRange[1] + float(k)/10 * (maxRange[1]-minRange[1])
		print "    coord=", (xcoord,ycoord), "result=", uwpytools.FieldVariable_InterpolateValueAt(VelocityMagnitudeField, (xcoord,ycoord) )
	print "------------------------------------"
print ""
print ""
print ""

# continue through till end 
uwpytools.RunMainLoop()


uwpytools.Finalise()
