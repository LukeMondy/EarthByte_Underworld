#!/bin/sh
cd /mnt/scratch/uw2_spherical/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/julian/scratch/uw2_spherical/libUnderworld/build/bin/StGermain LidDrivenFEM.xml /mnt/scratch/uw2_spherical/libUnderworld/Spherical/SysTests/Regression/expected/LidDrivenFEM-referenceTest/credo-analysis.xml 
