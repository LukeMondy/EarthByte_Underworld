#!/bin/sh
cd /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression
mpiexec -np 1 /home/jgiordani/codes/clean_uw/libUnderworld/build/bin/StGermain RTP3D_pic.xml /home/jgiordani/codes/clean_uw/libUnderworld/Spherical/SysTests/Regression/expected/RTP3D_pic-referenceTest/credo-analysis.xml 
