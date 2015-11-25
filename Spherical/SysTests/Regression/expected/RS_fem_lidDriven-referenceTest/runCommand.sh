#!/bin/sh
cd /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/jgiordani/codes/clean_uw/libUnderworld/build/bin/StGermain RS_fem_lidDriven.xml /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression/expected/RS_fem_lidDriven-referenceTest/credo-analysis.xml 
