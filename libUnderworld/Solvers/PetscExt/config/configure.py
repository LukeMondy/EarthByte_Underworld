#!/usr/bin/env python

import os
import sys
import time

import petscconf
import log
import check
import lapack

if not hasattr(sys, 'version_info') or not sys.version_info[1] >= 2 or not sys.version_info[0] >= 2:
  print '**** You must have Python version 2.2 or higher to run config/configure.py ******'
  print '*           Python is easy to install for end users or sys-admin.               *'
  print '*                   http://www.python.org/download/                             *'
  print '*                                                                               *'
  print '*            You CANNOT configure PETScExt without Python                       *'
  print '*********************************************************************************'
  sys.exit(4)
  

# support a few standard configure option types 
for l in range(1,len(sys.argv)):
  name = sys.argv[l]
  if name.startswith('--enable'):
    sys.argv[l] = name.replace('--enable','--with')
    if name.find('=') == -1: sys.argv[l] += '=1'
  if name.startswith('--disable'):
    sys.argv[l] = name.replace('--disable','--with')
    if name.find('=') == -1: sys.argv[l] += '=0'
    elif name.endswith('=1'): sys.argv[l].replace('=1','=0')
  if name.startswith('--without'):
    sys.argv[l] = name.replace('--without','--with')
    if name.find('=') == -1: sys.argv[l] += '=0'
    elif name.endswith('=1'): sys.argv[l].replace('=1','=0')


# Check if enviroment is ok
print 'Checking environment...'
if 'PETSCEXT_DIR' not in os.environ:
  sys.exit('ERROR: PETSCEXT_DIR enviroment variable is not set')
petscextdir = os.environ['PETSCEXT_DIR']
petscarch = os.environ['PETSC_ARCH']

# Check PETSCEXT_DIR mathches current working directory
cwd = os.getcwd()
if not cwd == petscextdir:
  error_mesg = os.sep.join( ['ERROR: PETSCEXT_DIR must match current working directory: Yours was', petscextdir,', should be: '])
  error_mesg = os.sep.join([error_mesg,cwd])
  sys.exit( error_mesg )

if not os.path.exists(petscextdir) or not os.path.exists(os.sep.join([petscextdir,'conf'])):
  sys.exit('ERROR: PETSCEXT_DIR enviroment variable is not valid as petscextdir/conf does not exist')
os.chdir(petscextdir)
if 'PETSC_DIR' not in os.environ:
  sys.exit('ERROR: PETSC_DIR enviroment variable is not set [2]')
petscdir = os.environ['PETSC_DIR']
#if not os.path.exists(petscdir) or not os.path.exists(os.sep.join([petscdir,petscarch])):
#  sys.exit('ERROR: PETSC_DIR enviroment variable is not valid [3]')

# Setup installation directory
prefixdir = ''
#  def configureInstall(self):
#if self.framework.argDB['prefix']:
#  prefixdir = self.framework.argDB['prefix']
#else:
# prefixdir = petscdir

for arg in sys.argv[1:]:
  L = arg.count ( '--prefix=' )
  if L == 1:
    prefixdir = arg.replace ( '--prefix=', '' )



# Check some information about PETSc configuration
petscconf.Load(petscdir)
if petscconf.VERSION < '3.1.0':
  sys.exit('ERROR: This PETScExt version is not compatible with PETSc version '+petscconf.VERSION) 
if not petscconf.PRECISION in ['double','single','matsingle']:
  sys.exit('ERROR: This PETScExt version does not work with '+petscconf.PRECISION+' precision')

# Create architecture directory and configuration file
archdir = os.sep.join([petscextdir,petscconf.ARCH])
if not os.path.exists(archdir):
  try:
    os.mkdir(archdir)
  except:
    sys.exit('ERROR: cannot create architecture directory ' + archdir)
try:
  petscextconf = open(os.sep.join([archdir,'petscextconf']),'w')
  if not prefixdir:
    prefixdir = petscextdir
  petscextconf.write('PETSCEXT_INSTALL_DIR =' + prefixdir +'\n')
except:
  sys.exit('ERROR: cannot create configuration file in ' + archdir)


# Create conf dir
confdir = os.sep.join([petscextdir,petscconf.ARCH,'conf'])
if not os.path.exists(confdir):
  try:
    os.mkdir(confdir)
  except:
    sys.exit('ERROR: cannot create conf directory ' + confdir)



# Open log file
#log.Open( 'configure-' + petscconf.ARCH + '.log' )
# Write configure.log to PETSCEXT_DIR/PETSC_ARCH/conf/configure.log
log.Open( confdir + '/configure.log' )
log.Write('='*80)
log.Write('Starting Configure Run at '+time.ctime(time.time()))
log.Write('Configure Options: '+str.join(' ',sys.argv))
log.Write('Working directory: '+os.getcwd())
log.Write('Python version:\n' + sys.version)
log.Write('='*80)

# Check if PETSc is working
log.Println('Checking PETSc installation...')
if petscconf.VERSION > '3.1.0':
  log.Println('WARNING: PETSc version '+petscconf.VERSION+' is newer than PETScExt version')
if petscconf.RELEASE != '1':
  log.Println('WARNING: using PETSc development version')
#if petscconf.INSTALL_DIR != petscdir:
#  log.Println('WARNING: PETSC_DIR does not point to PETSc installation path')
#  log.Println('  PETSC_DIR='+petscdir+'') 
#  log.Println('  INSTALL_DIR='+petscconf.INSTALL_DIR+'')
if not check.Link([],[],[]):
  log.Exit('ERROR: Unable to link with PETSc')


# Check for missing LAPACK functions
missing = lapack.Check(petscextconf)

petscextconf.close()

log.Println('')
log.Println('='*80)
log.Println('PETScExt Configuration')
log.Println('='*80)
log.Println('')
log.Println('PETScExt source directory:')
log.Println(' '+petscextdir)
log.Println('PETScExt install directory:')
log.Println(' '+prefixdir)  
log.Println('PETSc directory:')
log.Println(' '+petscdir)
log.Println('Architecture "'+petscconf.ARCH+'" with '+petscconf.PRECISION+' precision '+petscconf.SCALAR+' numbers')
if petscconf.MPIUNI:
  log.Println('  Uniprocessor version without MPI')
if missing:
  log.Println('LAPACK missing functions:')
  log.Print('  ')
  for i in missing: log.Print(i)
  log.Println('')
  log.Println('')
  log.Println('WARNING: Some PETScExt functionality will not be available')
  log.Println('PLEASE reconfigure and recompile PETSc with a full LAPACK implementation')
print



# Open petscext_variable file
# Write variables files to PETSCEXT_DIR/PETSC_ARCH/conf/petscext_variables
petscext_variables = open( confdir + '/petscext_variables', 'w' )
petscext_variables.write('PETSCEXT_HELPERS_LIB_BASIC = -lpetscext\n')
petscext_variables.write('PETSCEXT_VEC_LIB_BASIC = -lpetscext\n')
petscext_variables.write('PETSCEXT_MAT_LIB_BASIC = -lpetscext\n')
petscext_variables.write('PETSCEXT_KSP_LIB_BASIC = -lpetscext\n')
petscext_variables.write('PETSCEXT_SNES_LIB_BASIC = -lpetscext\n')
petscext_variables.write('PETSCEXT_UTILS_LIB_BASIC = -lpetscext\n')
petscext_variables.write('LIBNAME = ${INSTALL_LIB_DIR}/libpetscext.${AR_LIB_SUFFIX}\n')
petscext_variables.close()
