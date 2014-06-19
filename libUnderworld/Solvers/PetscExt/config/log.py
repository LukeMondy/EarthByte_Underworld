import sys

import petscconf

def Open(filename):
  global f
  f = open(filename,'w')
  return

def Println(string):
  print string
  f.write(string)
  f.write('\n')

def Print(string):
  print string,
  f.write(string+' ')
  
def Write(string):
  f.write(string)
  f.write('\n')
  
def Exit(string):
  f.write(string)
  f.write('\n')
  f.close()
  print string
  sys.exit('ERROR: See "configure-' + petscconf.ARCH + '.log' '" file for details')

