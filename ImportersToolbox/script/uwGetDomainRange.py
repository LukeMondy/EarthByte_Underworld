#!/usr/bin/env python

import os, sys, errno
import getopt
import string
import re

def fileType( filename ):

   fileopened = open( filename, 'r' )
   isGocad=False
   isGeomodeller=False
   isVTK=False
   gocad=0
   geomodeller=0
   vtk=0
   for line in fileopened:
      if ((re.search("^\s*AXIS_O\s", line) !=None) | (re.search("^\s*AXIS_U\s"  , line) !=None) | (re.search("^\s*AXIS_V\s"   , line) !=None) |
          (re.search("^\s*AXIS_W\s", line) !=None) | (re.search("^\s*AXIS_MIN\s", line) !=None) | (re.search("^\s*AXIS_MAX\s" , line) !=None) |
          (re.search("^\s*AXIS_N\s", line) !=None) | (re.search("^\s*AXIS_D\s"  , line) !=None) | (re.search("^\s*ZPOSITIVE\s", line) !=None) ):
         gocad+=1
      if ((re.search("^\s*nx\s", line) !=None) | (re.search("^\s*ny\s", line) !=None) | (re.search("^\s*nz\s", line) !=None) |
          (re.search("^\s*x0\s", line) !=None) | (re.search("^\s*y0\s", line) !=None) | (re.search("^\s*z0\s", line) !=None) |
          (re.search("^\s*dx\s", line) !=None) | (re.search("^\s*dy\s", line) !=None) | (re.search("^\s*dz\s", line) !=None) ) :
         geomodeller+=1
      if ((re.search("^\s*DATASET\s", line) !=None) | (re.search("^\s*DIMENSIONS\s", line) !=None) | 
          (re.search("^\s*ORIGIN\s", line) !=None)  | (re.search("^\s*SPACING\s", line) !=None)    | (re.search("^\s*POINT_DATA\s", line) !=None) ):
         vtk+=1
      # do this check so that we don't read the entire file
      if geomodeller > 7:
         break
      if vtk > 3:
         break
   fileopened.close()
   if gocad > 4:
      isGocad=True
   if geomodeller > 4:
      isGeomodeller=True
   if vtk > 3:
      isVTK=True

   if( isGocad & isGeomodeller & isVTK ):
      print "Error... Unable to determine if file is Gocad or Geomodeller or VTK."
      sys.exit( 2 )

   if isGocad:
      return "Gocad"
   if isGeomodeller:
      return "Geomodeller"
   if isVTK:
      return "VTK"



def getDomainGocad( filename, eps ):

   AXIS_O=None
   AXIS_U=None
   AXIS_V=None
   AXIS_W=None
   AXIS_MIN=None
   AXIS_MAX=None
   AXIS_N=None
   AXIS_D=None
   CELL_WIDTH=[0,0,0]
   
   fileopened = open( filename, 'r' )
   for line in fileopened:
      if   re.search("^\s*AXIS_O\s", line):
         AXIS_O    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_U\s", line):
         AXIS_U    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_V\s", line):
         AXIS_V    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_W\s", line):
         AXIS_W    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_MIN\s", line):
         AXIS_MIN  = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_MAX\s", line):
         AXIS_MAX  = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_N\s", line):
         AXIS_N    = [ int(string.split(line)[1]), int(string.split(line)[2]), int(string.split(line)[3]) ]
      elif re.search("^\s*AXIS_D\s", line):
         AXIS_D    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
      elif re.search("^\s*ZPOSITIVE\s", line):
         ZPOSITIVE = [ string.split(line)[1] ]
   fileopened.close()
   
   if( (AXIS_O==None) | (AXIS_U==None) | (AXIS_V==None) | (AXIS_W==None) | (AXIS_MIN==None) | (AXIS_MAX==None) ): 
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )
   
#   # first check for orthogonal & axis aligned axis
#   dot = 0
#   for item in AXIS_U:
#      dot += pow(item,2)
#   if ( pow(sum(AXIS_U),2) < 0.9999999*dot ):
#      print "Axis do not appear to be aligned with x,y,z axis. We cannot currently handle these voxets"
#      sys.exit( 2 )
#   dot = 0
#   for item in AXIS_V:
#      dot += pow(item,2)
#   if ( pow(sum(AXIS_V),2) < 0.9999999*dot ):
#      print "Axis do not appear to be aligned with x,y,z axis. We cannot currently handle these voxets"
#      sys.exit( 2 )
#   dot = 0
#   for item in AXIS_W:
#      dot += pow(item,2)
#   if ( pow(sum(AXIS_W),2) < 0.9999999*dot ):
#      print "Axis do not appear to be aligned with x,y,z axis. We cannot currently handle these voxets"
#      sys.exit( 2 )

   # first check for orthogonal & axis aligned axis
   if ( (AXIS_U[1]!=0) | (AXIS_U[2]!=0) ):
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )
   if ( (AXIS_V[0]!=0) | (AXIS_V[2]!=0) ):
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )
   if ( (AXIS_W[0]!=0) | (AXIS_W[1]!=0) ):
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )

   x1 = AXIS_O[0] + AXIS_U[0]*AXIS_MIN[0]
   x2 = AXIS_O[0] + AXIS_U[0]*AXIS_MAX[0]   
   y1 = AXIS_O[1] + AXIS_V[1]*AXIS_MIN[1]
   y2 = AXIS_O[1] + AXIS_V[1]*AXIS_MAX[1]   
   z1 = AXIS_O[2] + AXIS_W[2]*AXIS_MIN[2]
   z2 = AXIS_O[2] + AXIS_W[2]*AXIS_MAX[2]   

   if (ZPOSITIVE == "Depth"):
      z1 = -z1
      z2 = -z2
   minx = min(x1,x2)
   maxx = max(x1,x2)
   miny = min(y1,y2)
   maxy = max(y1,y2)
   minz = min(z1,z2)
   maxz = max(z1,z2)

   if AXIS_N:
      CELL_WIDTH[0] = (maxx-minx)/float(AXIS_N[0]-1)
      CELL_WIDTH[1] = (maxy-miny)/float(AXIS_N[1]-1)
      CELL_WIDTH[2] = (maxz-minz)/float(AXIS_N[2]-1)
   elif AXIS_D:
      CELL_WIDTH[0] = AXIS_D[0]
      CELL_WIDTH[1] = AXIS_D[1]
      CELL_WIDTH[2] = AXIS_D[2]
   else:
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )

   return [ minx-CELL_WIDTH[0]/2, maxx+CELL_WIDTH[0]/2, miny-CELL_WIDTH[1]/2, maxy+CELL_WIDTH[1]/2, minz-CELL_WIDTH[2]/2, maxz+CELL_WIDTH[2]/2 ]



def getDomainGeomodeller( filename, eps ):

   AXIS_O=[None,None,None]
   AXIS_U=[1,0,0]
   AXIS_V=[0,1,0]
   AXIS_W=[0,0,1]
   AXIS_MIN=[0,0,0]
   AXIS_MAX=[1,1,1]
   AXIS_N=[None,None,None]
   AXIS_D=[None,None,None]
   CELL_WIDTH=[None,None,None]
   
   count=0
   fileopened = open( filename, 'r' )
   for line in fileopened:
      if   re.search("^\s*x0\s", line):
         AXIS_O[0] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*y0\s", line):
         AXIS_O[1] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*z0\s", line):
         AXIS_O[2] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*nx\s", line):
         AXIS_N[0] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*ny\s", line):
         AXIS_N[1] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*nz\s", line):
         AXIS_N[2] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*dx\s", line):
         AXIS_D[0] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*dy\s", line):
         AXIS_D[1] = float(string.split(line)[1])
         count+=1
      if   re.search("^\s*dz\s", line):
         AXIS_D[2] = float(string.split(line)[1])
         count+=1
      if count == 9:
         break
   fileopened.close()

   x1 = AXIS_O[0]
   x2 = AXIS_O[0] + AXIS_D[0]*float((AXIS_N[0]-1))
   y1 = AXIS_O[1]
   y2 = AXIS_O[1] + AXIS_D[1]*float((AXIS_N[1]-1))
   z1 = AXIS_O[2]
   z2 = AXIS_O[2] + AXIS_D[2]*float((AXIS_N[2]-1))

   minx = min(x1,x2)
   maxx = max(x1,x2)
   miny = min(y1,y2)
   maxy = max(y1,y2)
   minz = min(z1,z2)
   maxz = max(z1,z2)

   CELL_WIDTH[0] = AXIS_D[0]
   CELL_WIDTH[1] = AXIS_D[1]
   CELL_WIDTH[2] = AXIS_D[2]

   return [ minx-CELL_WIDTH[0]/2, maxx+CELL_WIDTH[0]/2, miny-CELL_WIDTH[1]/2, maxy+CELL_WIDTH[1]/2, minz-CELL_WIDTH[2]/2, maxz+CELL_WIDTH[2]/2 ]

def getDomainVTK( filename, eps ):
   AXIS_O=[None,None,None]
   AXIS_U=[1,0,0]
   AXIS_V=[0,1,0]
   AXIS_W=[0,0,1]
   AXIS_MIN=[0,0,0]
   AXIS_MAX=[1,1,1]
   AXIS_N=[None,None,None]
   AXIS_D=[None,None,None]
   CELL_WIDTH=[None,None,None]
   
   fileopened = open( filename, 'r' )
   check = 0
   for line in fileopened:
      if   re.search("^\s*ORIGIN\s", line):
         AXIS_O    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
         check+=1
      elif re.search("^\s*DIMENSIONS\s", line):
         AXIS_N    = [ int(string.split(line)[1]), int(string.split(line)[2]), int(string.split(line)[3]) ]
         check+=1
      elif re.search("^\s*SPACING\s", line):
         AXIS_D    = [ float(string.split(line)[1]), float(string.split(line)[2]), float(string.split(line)[3]) ]
         check+=1
      if check == 3:
         break
   fileopened.close()
   
   if( (AXIS_O==None) | (AXIS_N==None) | (AXIS_D==None) ): 
      print "NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE NOTAVAILABLE"
      sys.exit( 2 )
   
   x1 = AXIS_O[0]
   x2 = AXIS_O[0] + AXIS_D[0]*float((AXIS_N[0]-1))
   y1 = AXIS_O[1]
   y2 = AXIS_O[1] + AXIS_D[1]*float((AXIS_N[1]-1))
   z1 = AXIS_O[2]
   z2 = AXIS_O[2] + AXIS_D[2]*float((AXIS_N[2]-1))

   minx = min(x1,x2)
   maxx = max(x1,x2)
   miny = min(y1,y2)
   maxy = max(y1,y2)
   minz = min(z1,z2)
   maxz = max(z1,z2)

   CELL_WIDTH[0] = AXIS_D[0]
   CELL_WIDTH[1] = AXIS_D[1]
   CELL_WIDTH[2] = AXIS_D[2]

   return [ minx-CELL_WIDTH[0]/2, maxx+CELL_WIDTH[0]/2, miny-CELL_WIDTH[1]/2, maxy+CELL_WIDTH[1]/2, minz-CELL_WIDTH[2]/2, maxz+CELL_WIDTH[2]/2 ]


def usage():
   #      12345678901234567890123456789012345678901234567890123456789012345678901234567890
   print "Tool to simply retrieve a gocad or geomodeller voxel domain range"
   print "Results are sent to standard out as \"minx maxx miny maxy minz maxz\" "
   print "This script will attempted to automatically determine if the file is Geomodeller or Gocad."
   print ""
   print "usage: uwGetGocadDomainRange.py --filename=FooGocadVoxet.vo"
   print ""
   print "Copyright (c) Monash University, 2012"
   print ""
   print "--help (-e)                            Prints out this help information."
   print "--filename=name (-n name)              Specify the gocad voxet filename."
   print "--domFactor=    (-d name)              Specify a factor to increase (domain in all directions) by."
   print "--format=       (-f name)              Specify 'Gocad' or 'Geomodeller' format."
   print ""


def main():
   # Parse options...
   try:
      opts, args = getopt.getopt(sys.argv[1:], "hn:d:f:", ["help", "filename=", "domFactor=", "format="])
   except getopt.GetoptError, err:
      # print help information and exit:
      print str( err ) # will print something like "option -a not recognized"
      usage()
      sys.exit( 2 )
   name = None
   eps = 0.0
   format = None
   for o, a in opts:
      if o in ("-h", "--help"):
         usage()
         sys.exit()
      elif o in ("-n", "--filename"):
         name = a
      elif o in ("-d", "--domFactor"):
         eps = float(a)
      elif o in ("-f", "--format"):
         format = a
      else:
         assert False, "unhandled option"

   if( (name == None) | (name =="")):
      print "VoxetMinX VoxetMaxX VoxetMinY VoxetMaxY VoxetMinZ VoxetMaxZ"
      sys.exit( 2 )

   if( not os.path.isfile(name) ):
      print "File not found:",name
      sys.exit( 2 )

   if format == None:
      format = fileType( name )
   if(format=="Gocad"):
      dom = getDomainGocad( name, eps )
   if(format=="Geomodeller"):
      dom = getDomainGeomodeller( name, eps )
   if(format=="VTK"):
      dom = getDomainVTK( name, eps )

   sizeI = dom[1]-dom[0]
   sizeJ = dom[3]-dom[2]
   sizeK = dom[5]-dom[4]

   print dom[0] - eps*sizeI, dom[1] + eps*sizeI, dom[2] - eps*sizeJ, dom[3] + eps*sizeJ, dom[4] - eps*sizeK, dom[5] + eps*sizeK


if __name__ == "__main__":
   main()

