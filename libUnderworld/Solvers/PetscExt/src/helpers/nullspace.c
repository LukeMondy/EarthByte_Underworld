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

#include "petscext_helpers.h"


/*

A set of simple functions to help you 
    i) estimate the size of the constant null space
   ii) detect the precesence of the null space using a default value (PETSCEXT_NULLSPACE_SIZE)
  iii) strip the constant null space from a vector

*/


#undef __FUNCT__
#define __FUNCT__ "PetscExtMatComputeConstantNullSpace"
PetscErrorCode PetscExtMatComputeConstantNullSpace( Mat A, PetscReal *nrm )
{
	PetscInt N;
	PetscScalar sum;
	Vec l,r;
	
	
	PetscFunctionBegin;
	
	MatGetVecs( A, &r, &l ); // l = A r
	
	
	VecGetSize(r,&N);
	sum  = 1.0/N;
	VecSet(r,sum);
	
	/* {l} = [A] {r} */
	MatMult( A,r, l );
	
	VecNorm(l,NORM_2,nrm);
	
	Stg_VecDestroy( &l );
	Stg_VecDestroy( &r );
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscExtMatContainsConstantNullSpace"
PetscErrorCode PetscExtMatContainsConstantNullSpace( Mat A, PetscTruth *has_cnst_nullsp )
{
	PetscInt N;
	PetscScalar sum;
	PetscReal nrm;
	Vec l,r;
	
	
	PetscFunctionBegin;
	
	MatGetVecs( A, &r, &l ); // l = A r
	
	
	VecGetSize(r,&N);
	sum  = 1.0/N;
	VecSet(r,sum);
	
	/* {l} = [A] {r} */
	MatMult( A,r, l );
	
	VecNorm(l,NORM_2,&nrm);
	if (nrm < PETSCEXT_NULLSPACE_SIZE) {
		*has_cnst_nullsp = PETSC_TRUE;
	}
	else {
		*has_cnst_nullsp = PETSC_FALSE;
	}
	
	Stg_VecDestroy( &l );
	Stg_VecDestroy( &r );
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PetscExtVecRemoveConstantNullspace"
PetscErrorCode PetscExtVecRemoveConstantNullspace( Vec v )
{
	PetscInt N;
	PetscScalar sum;
	
	
	PetscFunctionBegin;
	
	VecGetSize( v, &N );
	if( N > 0 ) {
		VecSum( v, &sum );
		sum  = sum/( -1.0*N );
		VecShift( v, sum );
	}
	
	PetscFunctionReturn(0);
}



