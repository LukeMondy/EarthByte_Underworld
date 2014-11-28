#!/usr/bin/env python
'''
This example shows how you might run multiple related (or not) simulations from a single script. 
Here, StGermain is fully constructed and torn down between simulations.  Note that only one
call is made to MPI_Init, so this is useful for running on large hpc facilities under a single job launch. 

Here we simply modify the density between simulations.  This is effected by changing the input dictionary. 
'''


import underworld
import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 1:
	print os.environ['UWPATH']

print "======================="

for ii in range(1,5):
	# init using json
	underworld.Init("RayleighTaylorBenchmark.json")

	# grab the dict
	stgdict = underworld.GetCurrentDictionary()

	# create a distinct path for each sim
	stgdict["outputPath"] = "outputPath" + "_" + str(ii)

	# modify the density for each sim
	stgdict["components"]["lightLayer"]["density"] = float(ii)/10

	# don't forget to set the dict back again to effect the above changes
	underworld.SetDictionary(stgdict)

	underworld.Construct()
	# run through till end 
	underworld.RunMainLoop()
	underworld.Finalise()
