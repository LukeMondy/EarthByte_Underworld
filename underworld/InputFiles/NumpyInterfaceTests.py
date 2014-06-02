# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import uwpytools
import numpy as np
import math

# <codecell>

uwpytools.Init()

# <codecell>

# init using json
uwpytools.InitWithArgs("RayleighTaylorBenchmark.json")

# <codecell>

uwpytools.Construct()
uwpytools.BuildAndInitialise()

# <codecell>

velthing = uwpytools.GetLiveComponent("VelocityField")

# <codecell>

mesh = uwpytools.GetLiveComponent("linearMesh")

# <markdowncell>

# # lets test out the operator overloading

# <codecell>

velthing[0]

# <codecell>

mesh[0:4225:1000]

# <markdowncell>

# # we can also set values using the numpy interface

# <codecell>

velthing[10][0] = 1

# <codecell>

velthing[10]

# <codecell>

mesh[0][0] = 0.0001

# <codecell>

mesh[0:4225:1000]

# <markdowncell>

# # great success..  now, lets do some further testing.  first, lets check the original implementation benchmark:

# <codecell>

# the old interface version.  note that this function must write into stg memory
def TrigICOrig(meshGuy, velGuy):
    import uwpytools.c_arrays as c_arrays
    valuePtr = c_arrays.DoubleArray(2)

    for ii in range(0,uwpytools.StgDomain.Mesh_GetDomainSize(meshGuy,0)):
    	position = c_arrays.DoubleArray_frompointer(uwpytools.StgDomain.Mesh_GetVertex( meshGuy, ii ))
        valuePtr[0] = math.sin(position[0])**2
        valuePtr[1] = math.cos(position[1])**2
        uwpytools.StgFEM.FeVariable_SetValueAtNode( velGuy, ii, valuePtr.cast() )

# <codecell>

timeit TrigICOrig(mesh, velthing)

# <markdowncell>

# # so, 34ms.  let's try using the operator overload method

# <codecell>

# the overloaded square bracket interface function
def TrigICOverload(meshGuy, velGuy):
    for ii in range(0, velGuy.getAsNumpyArray().shape[0]):
        velGuy[ii][0] = math.sin(meshGuy[ii][0])**2
        velGuy[ii][1] = math.cos(meshGuy[ii][1])**2

# <codecell>

# grab a copy to compare with later
velArrayCopy = np.copy(velthing.getAsNumpyArray())

# <codecell>

timeit TrigICOverload(mesh, velthing)

# <codecell>

# check they are the same
np.linalg.norm(velArrayCopy-velthing.getAsNumpyArray()) 

# <markdowncell>

# # so, only minor improvements in efficiency here at 29ms.. was expecting more.  in any case, certainly huge improvements in user interface. let's test out directly manipulating numpy arrays.  first, retrieve arrays from stg components:

# <codecell>

velArray = velthing.getAsNumpyArray()

# <codecell>

meshArray = mesh.getAsNumpyArray()

# <markdowncell>

# # great, now lets define a function to directly utilse these arrays:

# <codecell>

# direct numpy function
def TrigICNumpy1(meshArray, fieldArray):
    for meshRow, fieldRow in zip(meshArray, fieldArray):
        fieldRow[0] = math.sin(meshRow[0])**2
        fieldRow[1] = math.cos(meshRow[1])**2

# <markdowncell>

# # lets test this function out..  note that it directly writes to the stg memory via the numpy proxy

# <codecell>

timeit TrigICNumpy1(meshArray, velArray)

# <codecell>

# check they are the same
np.linalg.norm(velArrayCopy-velArray) 

# <markdowncell>

# # ok, factor of 3-4 faster.  not bad.  how else can we write this function?

# <codecell>

# this version does the same but is faster
def TrigICNumpy2(meshArray, fieldArray):
    fieldArray[:,0] = np.sin(meshArray[:,0])**2
    fieldArray[:,1] = np.cos(meshArray[:,1])**2

# <codecell>

timeit TrigICNumpy2(meshArray, velArray)

# <codecell>

# check they are the same
np.linalg.norm(velArrayCopy-velArray) 

# <markdowncell>

# # a further factor of 60 improvement. great success.. at least for simplish functions.  lets now try using the external package numexpr, which should provide for similar improvements even for complex functions

# <codecell>

# ok, lets try doing this using numexpr
def TrigICNumexpr(meshArray,velArray):
    import numexpr as ne
    #  numexpr doesnt seem to handle slices
    meshslicex = meshArray[:,0]
    meshslicey = meshArray[:,1]
    ne.evaluate( "sin(meshslicex)**2", out=velArray[:,0])
    ne.evaluate( "cos(meshslicey)**2", out=velArray[:,1])

# <codecell>

timeit TrigICNumexpr(meshArray, velArrayCopy)

# <codecell>

# check they are the same... note that velArray has been overwritten using the above, so compare to copy
np.linalg.norm(velArrayCopy-velArray)

# <markdowncell>

# # ok, this one is slower than than regular numpy version.  remains to see how it does where the complexity is stepped up.

# <codecell>

uwpytools.Finalise()

