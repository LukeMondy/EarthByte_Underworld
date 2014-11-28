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

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "petscext_helpers.h"

#undef __FUNCT__  
#define __FUNCT__ "MatNullSpaceQuery"
PetscErrorCode PETSCMAT_DLLEXPORT MatNullSpaceQuery(Mat mat,PetscReal *cnst_nullsp_val,PetscTruth *has_cnst_nullsp,PetscReal nullsp_val[],PetscTruth has_nullsp[])
{
  PetscScalar    sum;
  PetscReal      nrm;
  PetscInt       j,n,N,m;
  PetscErrorCode ierr;
  Vec            l,r;
  Vec            vec;
  PetscViewer    viewer;
  MatNullSpace   sp;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mat,MAT_COOKIE,2);
  sp = mat->nullsp;
  PetscValidHeaderSpecific(sp,MAT_NULLSPACE_COOKIE,1);

  if (!sp) {
    Stg_SETERRQ(PETSC_ERR_ORDER,"Mat has no null space attached. Must call MatNullSpaceAttach() first.");
  }

  n = sp->n;
#if ( (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5) )
  if (n) {
    ierr = VecDuplicate(sp->vecs[0],&vec);CHKERRQ(ierr);
  } else {
    ierr = MatGetLocalSize(mat,&m,PETSC_NULL);CHKERRQ(ierr);
    ierr = VecCreateMPI(((PetscObject)sp)->comm,m,PETSC_DETERMINE,&vec);CHKERRQ(ierr);
  }
  l = vec;
#else
  if (!sp->vec) {
    if (n) {
      ierr = VecDuplicate(sp->vecs[0],&sp->vec);CHKERRQ(ierr);
    } else {
      ierr = MatGetLocalSize(mat,&m,PETSC_NULL);CHKERRQ(ierr);
      ierr = VecCreateMPI(((PetscObject)sp)->comm,m,PETSC_DETERMINE,&sp->vec);CHKERRQ(ierr);
    }
  }
  l = sp->vec;
#endif

  ierr = PetscViewerASCIIGetStdout(((PetscObject)sp)->comm,&viewer);CHKERRQ(ierr);

  if (sp->has_cnst) {
    ierr = VecDuplicate(l,&r);CHKERRQ(ierr);
    ierr = VecGetSize(l,&N);CHKERRQ(ierr);
    sum  = 1.0/N;
    ierr = VecSet(l,sum);CHKERRQ(ierr);
    ierr = MatMult(mat,l,r);CHKERRQ(ierr);
    ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
    if (nrm < 1.e-7) {
      if (has_cnst_nullsp){ *has_cnst_nullsp = PETSC_TRUE; }
    }
    else {
      if (has_cnst_nullsp){ *has_cnst_nullsp = PETSC_FALSE; }
    }
    if (cnst_nullsp_val) { *cnst_nullsp_val = nrm; }
    ierr = Stg_VecDestroy( &r);CHKERRQ(ierr);
  }

  for (j=0; j<n; j++) {
    ierr = (*mat->ops->mult)(mat,sp->vecs[j],l);CHKERRQ(ierr);
    ierr = VecNorm(l,NORM_2,&nrm);CHKERRQ(ierr);
    if (nrm < 1.e-7) {
      if (has_nullsp){ has_nullsp[j] = PETSC_TRUE; }
    }
    else {
      if (has_nullsp){ has_nullsp[j] = PETSC_FALSE; }
    }
    if (nullsp_val) { nullsp_val[j] = nrm; }
  }

  if (sp->remove){
    Stg_SETERRQ(PETSC_ERR_SUP,"Cannot test a null space provided as a function with MatNullSpaceSetFunction()");
  }
  
  PetscFunctionReturn(0);
}
