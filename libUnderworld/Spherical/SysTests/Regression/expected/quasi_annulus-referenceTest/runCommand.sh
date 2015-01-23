#!/bin/sh
cd /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/jgiordani/codes/clean_uw/libUnderworld/build/bin/StGermain quasi_annulus.xml /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression/expected/quasi_annulus-referenceTest/credo-analysis.xml 
