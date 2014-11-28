import os
import sys

import petscconf
import log
import check

def Check(conf):
  log.Write('='*80)
  log.Println('Checking LAPACK library...')

  # =======================================================
  #  I've added wrappers for the functions listed below 
  # =======================================================
  functions = ['ppsv']

  if petscconf.SCALAR == 'real':
    if petscconf.PRECISION == 'single':
      prefix = 's'
    else:
      prefix = 'd'
  else:
    if petscconf.PRECISION == 'single':
      prefix = 'c'
    else:
      prefix = 'z'
  
  missing = []
  conf.write('PETSCEXT_MISSING_LAPACK =')
  
  for i in functions:
    f =  '#if defined(PETSC_HAVE_FORTRAN_UNDERSCORE) || defined(PETSC_BLASLAPACK_UNDERSCORE)\n'
    f += prefix + i + '_\n'
    f += '#elif defined(PETSC_HAVE_FORTRAN_CAPS)\n'
    f += prefix.upper() + i.upper() + '\n'
    f += '#else\n'
    f += prefix + i + '\n'
    f += '#endif\n'
   
    log.Write('=== Checking LAPACK '+prefix+i+' function...')
    if not check.Link([f],[],[]):
      missing.append(prefix + i)
      conf.write(' -DPETSCEXT_MISSING_LAPACK_' + i.upper())

    log.Write('=== Checking LAPACK '+i+' function...')
    if not check.Link([f],[],[]):
      missing.append(i)
      conf.write(' -DPETSCEXT_MISSING_LAPACK_' + i.upper())
  
  conf.write('\n')
  return missing
