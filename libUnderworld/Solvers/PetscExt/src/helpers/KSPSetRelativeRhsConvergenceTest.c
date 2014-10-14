/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2007 - 2008, 
 **        Dave A. May [dave.mayhem23@gmail.com]
 **        School of Mathematical Sciences, 
 **        Monash University
 **        Victoria, 3800
 **        Australia
 **
 **    Copyright (c) 2009 - 2010, 
 **        Dave A. May [dave.mayhem23@gmail.com]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Switzerland
 **
 **    Project:       PetscExt_v3.1
 **    Description:   A collection of extensions to the PETSc library.
 **
 **    All rights reserved.
 **    Redistribution and use in source and binary forms, with or without modification,
 **    are permitted provided that the following conditions are met:
 **
 **        * Redistributions of source code must retain the above copyright notice,
 **              this list of conditions and the following disclaimer.
 **        * Redistributions in binary form must reproduce the above copyright
 **              notice, this list of conditions and the following disclaimer in the
 **              documentation and/or other materials provided with the distribution.
 **        * Neither the name of the Monash University nor the names of its contributors
 **              may be used to endorse or promote products derived from this software
 **              without specific prior written permission.
 **
 **    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 **    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 **    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 **    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 **    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 **    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 **    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 **    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 **    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 **    OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#include <petscksp.h>
#include "private/compat/petsccompat.h"

#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3 ) )
  #include "petsc-private/kspimpl.h"   /*I "petscksp.h" I*/
#else
  #include "private/kspimpl.h"   /*I "petscksp.h" I*/
#endif


/*
Convergence occurs when
  rnorm < MAX ( rtol * norm_b, abstol )
Divergence detected when
  rnorm > dtol * rnorm_0

*/
#undef __FUNCT__
#define __FUNCT__ "KSPRelativeRhsConverged"
PetscErrorCode KSPRelativeRhsConverged(KSP ksp,PetscInt n,PetscReal rnorm,KSPConvergedReason *reason,void *dummy)
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	PetscValidHeaderSpecific(ksp,KSP_COOKIE,1);
	PetscValidPointer(reason,4);
	
	*reason = KSP_CONVERGED_ITERATING;
	
	if (!n) {
		PetscReal      snorm;
		
		if (ksp->normtype == KSP_NORM_UNPRECONDITIONED) {
			ierr = VecNorm(ksp->vec_rhs,NORM_2,&snorm);CHKERRQ(ierr);        /*     <- b'*b */
		}
		else if (ksp->normtype == KSP_NORM_PRECONDITIONED) {
			if( ksp->pc_side == PC_RIGHT ) {
				ierr = VecNorm(ksp->vec_rhs,NORM_2,&snorm);CHKERRQ(ierr);        /*     <- b'*b */
			}
			else {
				Vec z;
				ierr = VecDuplicate(ksp->vec_rhs,&z);CHKERRQ(ierr);
				ierr = KSP_PCApply(ksp,ksp->vec_rhs,z);CHKERRQ(ierr);
				ierr = VecNorm(z,NORM_2,&snorm);CHKERRQ(ierr);                 /*    dp <- b'*B'*B*b */
				ierr = Stg_VecDestroy( &z);CHKERRQ(ierr);
			}
		}
		else {
			Stg_SETERRQ(PETSC_ERR_SUP,"KSPRelativeRhsConvergenced() only supports KSPNormType = {KSP_NORM_UNPRECONDITIONED, KSP_NORM_PRECONDITIONED} ");
		}
		ksp->rnorm0 = snorm;
		ksp->ttol   = PetscMax(ksp->rtol*ksp->rnorm0,ksp->abstol);
	}
	
	/* Standard petsc test from KSPDefaultConverged() */
	if (rnorm != rnorm) {
		ierr = PetscInfo(ksp,"Linear solver has created a not a number (NaN) as the residual norm, declaring divergence \n");CHKERRQ(ierr);
		*reason = KSP_DIVERGED_NAN;
	} else if (rnorm <= ksp->ttol) {
		if (rnorm < ksp->abstol) {
			ierr = PetscInfo3(ksp,"Linear solver has converged. Residual norm %G is less than absolute tolerance %G at iteration %D\n",rnorm,ksp->abstol,n);CHKERRQ(ierr);
			*reason = KSP_CONVERGED_ATOL;
		} else {
			if (ksp->rnorm0) {
				ierr = PetscInfo4(ksp,"Linear solver has converged. Residual norm %G is less than relative tolerance %G times initial right hand side norm %G at iteration %D\n",rnorm,ksp->rtol,ksp->rnorm0,n);CHKERRQ(ierr);
			} else {
				ierr = PetscInfo4(ksp,"Linear solver has converged. Residual norm %G is less than relative tolerance %G times initial right hand side norm %G at iteration %D\n",rnorm,ksp->rtol,ksp->rnorm0,n);CHKERRQ(ierr);
			}
			*reason = KSP_CONVERGED_RTOL;
		}
	} else if (rnorm >= ksp->divtol*ksp->rnorm0) {
		ierr = PetscInfo3(ksp,"Linear solver is diverging. Initial right hand size norm %G, current residual norm %G at iteration %D\n",ksp->rnorm0,rnorm,n);CHKERRQ(ierr);
		*reason = KSP_DIVERGED_DTOL;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetRelativeRhsConvergenceTest"
PetscErrorCode KSPSetRelativeRhsConvergenceTest( KSP ksp )
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = KSPSetConvergenceTest( ksp, KSPRelativeRhsConverged, PETSC_NULL, PETSC_NULL );CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

