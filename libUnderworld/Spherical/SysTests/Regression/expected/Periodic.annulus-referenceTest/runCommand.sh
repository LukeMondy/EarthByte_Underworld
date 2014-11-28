#!/bin/sh
cd /mnt/scratch/clean_uw/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/julian/scratch/clean_uw/libUnderworld/build/bin/StGermain Periodic.annulus.xml /mnt/scratch/clean_uw/libUnderworld/Spherical/SysTests/Regression/expected/Periodic.annulus-referenceTest/credo-analysis.xml 
