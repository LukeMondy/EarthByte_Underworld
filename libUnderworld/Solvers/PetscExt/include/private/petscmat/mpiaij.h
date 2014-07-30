#if !defined(__MPIAIJ_H)
#define __MPIAIJ_H

#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==5) )
#include "mpiaij-3.5.0.h"
#endif
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==4) )
#include "mpiaij-3.4.5.h"
#endif
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==3) )
#include "mpiaij-3.3-p7.h"
#endif
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==2) )
#include "mpiaij-3.2-p7.h"
#endif
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==1) )
#include "mpiaij-3.1-p8.h"
#endif
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==0) )
#include "mpiaij-3.0.0-p12.h"
#endif

#endif
