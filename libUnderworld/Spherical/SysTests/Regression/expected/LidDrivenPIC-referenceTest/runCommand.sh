#!/bin/sh
cd /mnt/scratch/uw2_spherical/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/julian/scratch/uw2_spherical/libUnderworld/build/bin/StGermain LidDrivenPIC.xml /mnt/scratch/uw2_spherical/libUnderworld/Spherical/SysTests/Regression/expected/LidDrivenPIC-referenceTest/credo-analysis.xml 
