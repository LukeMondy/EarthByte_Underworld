#!/bin/sh
cd /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/jgiordani/codes/clean_uw/libUnderworld/build/bin/StGermain RS_pic.xml /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression/expected/RS_pic-referenceTest/credo-analysis.xml 
