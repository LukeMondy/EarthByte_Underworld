#!/usr/bin/env python

'''
This example shows how you might control the simulation flow, obtaining various metrics at certain timesteps.  
You might also modify various aspects of the simulation at different times in the simulation. 

'''

import underworld

# init using json
underworld.Init("RayleighTaylorBenchmark.xml")

underworld.Construct()

# grab the context
context = underworld._stgermain.GetLiveComponent("context")
# grab the velMag field
VelocityMagnitudeField = underworld._stgermain.GetLiveComponent("VelocityMagnitudeField")

minRange,maxRange = underworld.fields.tools.FieldVariable_GetMinAndMaxLocalCoords(VelocityMagnitudeField)

while context.timeStep < 2:
	underworld.Step()


print ""
print ""
print ""
if underworld.rank() == 0:
	print "------------------------------------"
	print "VelocityMagnitudeField local domain (min,max)", underworld.fields.tools.FieldVariable_GetMinAndMaxLocalCoords(VelocityMagnitudeField)
	print "VelocityMagnitudeField global domain (min,max)", underworld.fields.tools.FieldVariable_GetMinAndMaxGlobalCoords(VelocityMagnitudeField)
	print "------------------------------------"
print ""
print ""
# do some things
if underworld.rank() == 0:
	print "------------------------------------"
	print "At timestep", context.timeStep
	print "VelocityMagnitudeField min value", underworld.fields.tools.FieldVariable_GetMinFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField max value", underworld.fields.tools.FieldVariable_GetMaxFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField Integral", underworld.fields.tools.FeVariable_Integrate(VelocityMagnitudeField)

	print "VelocityMagnitudeField interpolation results:"
	for k in range(0,10):
		xcoord = minRange[0] + float(k)/10 * (maxRange[0]-minRange[0])
		ycoord = minRange[1] + float(k)/10 * (maxRange[1]-minRange[1])
		print "    coord=", (xcoord,ycoord), "result=", underworld.fields.tools.FieldVariable_InterpolateValueAt(VelocityMagnitudeField, (xcoord,ycoord) )
	print "------------------------------------"
print ""
print ""
print ""

# run for a few more steps
while context.timeStep < 4:
	underworld.Step()

print ""
print ""
print ""
if underworld.rank() == 0:
	print "------------------------------------"
	print "At timestep", context.timeStep
	print "VelocityMagnitudeField min value", underworld.fields.tools.FieldVariable_GetMinFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField max value", underworld.fields.tools.FieldVariable_GetMaxFieldMagnitude(VelocityMagnitudeField)
	print "VelocityMagnitudeField Integral", underworld.fields.tools.FeVariable_Integrate(VelocityMagnitudeField)

	print "VelocityMagnitudeField interpolation results:"
	for k in range(0,10):
		xcoord = minRange[0] + float(k)/10 * (maxRange[0]-minRange[0])
		ycoord = minRange[1] + float(k)/10 * (maxRange[1]-minRange[1])
		print "    coord=", (xcoord,ycoord), "result=", underworld.fields.tools.FieldVariable_InterpolateValueAt(VelocityMagnitudeField, (xcoord,ycoord) )
	print "------------------------------------"
print ""
print ""
print ""

# continue through till end 
underworld.RunMainLoop()


underworld.Finalise()
