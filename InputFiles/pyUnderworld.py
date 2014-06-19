#!/usr/bin/env python

''' 
This example can be used as a replacement for the regular StGermain/Underworld executable.  
You can feed your XML/JSON files directly in, along with the usual commandline additions. 
'''


import underworld

underworld.Init()
underworld.Construct()
underworld.RunMainLoop()
underworld.Finalise()
