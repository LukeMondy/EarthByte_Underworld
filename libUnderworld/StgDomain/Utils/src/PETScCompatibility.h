/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_PETScCompatibility_h__
#define __StgDomain_Utils_PETScCompatibility_h__

#include "petsc.h"

#if (((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=4)) || (PETSC_VERSION_MAJOR>3) )
   #define Stg_PCMGDefaultResidual NULL
#else
   #define Stg_PCMGDefaultResidual PCMGDefaultResidual
#endif

#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2 ) )
  	#define Stg_SETERRQ( arg1, arg2 ) SETERRQ( PETSC_COMM_SELF, arg1, arg2 )
  	#define Stg_SETERRQ1( arg1, arg2, arg3 ) SETERRQ2( PETSC_COMM_SELF, arg1, arg2, arg3 )
  	#define Stg_SETERRQ2( arg1, arg2, arg3, arg4 ) SETERRQ2( PETSC_COMM_SELF, arg1, arg2, arg3, arg4 )
  	#define Stg_MatZeroRows( arg1, arg2, arg3, arg4 ) MatZeroRows( arg1, arg2, arg3, arg4, NULL, NULL )
  	#define ISCreateGeneralWithArray( arg1, arg2, arg3, arg4 ) ISCreateGeneral( arg1, arg2, arg3, PETSC_USE_POINTER, arg4 )
  	#define PetscTruth PetscBool
  	#define PetscOptionsGetTruth PetscOptionsGetBool
  	#define PETSC_VIEWER_BINARY PETSCVIEWERBINARY
   #if ( PETSC_VERSION_MINOR > 2  )
      #define Stg_MatCreateAIJ MatCreateAIJ
      #define Stg_PetscObjectTypeCompare PetscObjectTypeCompare
   #else
      #define Stg_MatCreateAIJ MatCreateMPIAIJ
      #define Stg_PetscObjectTypeCompare PetscTypeCompare
   #endif
#else
/* need these for Mat Vec KSP etc in versions of petsc < 3.2 */
   #include "petscmat.h"
   #include "petscksp.h"
   #include "petscvec.h"
   #include "petscpc.h"
   #include "petscsnes.h"
   #define Stg_SETERRQ SETERRQ
   #define Stg_SETERRQ1 SETERRQ1
   #define Stg_SETERRQ2 SETERRQ2
   #define Stg_MatZeroRows MatZeroRows
   #define Stg_MatCreateAIJ MatCreateMPIAIJ
   #define Stg_PetscObjectTypeCompare PetscTypeCompare
#endif
/* wrapper functions for compatibility between Petsc version 3.2 and lower versions */
PetscErrorCode Stg_MatDestroy(Mat *A);
PetscErrorCode Stg_VecDestroy(Vec *A);
PetscErrorCode Stg_KSPDestroy(KSP *A);
PetscErrorCode Stg_PCDestroy(PC *A);
PetscErrorCode Stg_SNESDestroy(SNES *A);
PetscErrorCode Stg_VecScatterDestroy(VecScatter *ctx);
PetscErrorCode Stg_ISDestroy(IS *is);
PetscErrorCode Stg_PetscViewerDestroy(PetscViewer *viewer);
/** lets make Stg_MatLoad look the same as pre-petsc_version 3.2 MatLoad */
PetscErrorCode Stg_MatLoad(PetscViewer viewer, const MatType outtype,Mat *newmat);
PetscErrorCode Stg_VecLoad(PetscViewer viewer, const VecType outtype,Vec *newvec);


#endif /* __StgDomain_Utils_PETScCompatibility_h__ */

