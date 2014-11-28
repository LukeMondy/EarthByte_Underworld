
import os
import sys

def Load(petscdir):
  global ARCH,DIR,MAKE,SCALAR,PRECISION,MPIUNI,VERSION,RELEASE,INSTALL_DIR
  
  if 'PETSC_ARCH' in os.environ:
    ARCH = os.environ['PETSC_ARCH']
  else:
    try:
      f = open(os.sep.join([petscdir,ARCH,'conf','petscvariables']))
      ARCH = ''
      for l in f.readlines():
	if l.startswith('PETSC_ARCH_NAME='):
	  ARCH = l.split('=')[1].rstrip()
	  f.close()
	  break
      f.close()
    except:
      sys.exit('ERROR: PETSc must be configured first')
    if not ARCH:
      sys.exit('ERROR: please set enviroment variable PETSC_ARCH')

  MPIUNI = 0
  
  try:
    f = open(os.sep.join([petscdir,ARCH,'conf/petscvariables']))
    for l in f.readlines():
      (k,v) = l.split('=',1)
      k = k.strip()
      v = v.strip()
      if k == 'PETSC_SCALAR':
	SCALAR = v
      elif k == 'PETSC_PRECISION':
        PRECISION = v
      elif k == 'MPI_INCLUDE' and v.endswith('mpiuni'):
        MPIUNI = 1
      elif k == 'MAKE':
	MAKE = v
      elif k == 'INSTALL_DIR':
        INSTALL_DIR = v
    f.close()
  except:
#    sys.exit('ERROR: PETSc is not configured for architecture ' + ARCH)

        try:
                f = open(os.sep.join([petscdir,'conf/petscvariables']))
                for l in f.readlines():
                        (k,v) = l.split('=',1)
                        k = k.strip()
                        v = v.strip()
                        if k == 'PETSC_SCALAR':
                                SCALAR = v
                        elif k == 'PETSC_PRECISION':
                                PRECISION = v
                        elif k == 'MPI_INCLUDE' and v.endswith('mpiuni'):
                                MPIUNI = 1
                        elif k == 'MAKE':
                                MAKE = v
                        elif k == 'INSTALL_DIR':
                                INSTALL_DIR = v
                f.close()
        except:
                sys.exit('ERROR: PETSc is not configured for architecture ' + ARCH)


  try:
    f = open(os.sep.join([petscdir,'include','petscversion.h']))
    for l in f.readlines():
      l = l.split()
      if len(l) == 3:
        if l[1] == 'PETSC_VERSION_RELEASE':
	  RELEASE = l[2]
	if l[1] == 'PETSC_VERSION_MAJOR':
          major = l[2]
	elif l[1] == 'PETSC_VERSION_MINOR':
          minor = l[2]
	elif l[1] == 'PETSC_VERSION_SUBMINOR':
          subminor = l[2]
    f.close()
    VERSION = major + '.' + minor + '.' + subminor
  except:
    sys.exit('ERROR: file error while reading PETSC version')
