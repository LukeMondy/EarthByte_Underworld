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
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "petscext_helpers.h"

#undef __FUNCT__
#define __FUNCT__ "MatFreeComputeExplicitOperator"
PetscErrorCode MatFreeComputeExplicitOperator( Mat mf, MatType type, PetscReal tol, Mat *_A)
{
	PetscErrorCode ierr;
	PetscInt M,N,m,n;
	PetscInt s,e,c,i,*ridx;
	Mat A;
	Vec l,L,I_r;
	PetscScalar *vals, *_I_r;
	PetscMPIInt nnz_local,nnz,mierr;
	

	MatGetLocalSize( mf, &m,&n );
	MatGetSize( mf, &M,&N );
	
	MatCreate( ((PetscObject)mf)->comm, &A );
	MatSetSizes( A, m,n,M,N );
	MatSetType( A, type );
//	MatSetOption( A, MAT_COLUMNS_SORTED );
//	MatSetOption( A, MAT_ROWS_SORTED );
	
	MatGetVecs( mf, &I_r, &l );	/* {l} = mf {r} */
	VecDuplicate( l, &L );	

	MatGetOwnershipRange( A, &s,&e );
	PetscMalloc( sizeof(PetscScalar)*m, &ridx );
	for( i=0; i<m; i++ ) {
		ridx[i] = s+i;
	}
	
	nnz_local = 0;	
	for( c=0; c<N; c++ ) {
		
		VecZeroEntries( I_r );
		VecSetValue( I_r, c, 1.0, INSERT_VALUES );
		VecAssemblyBegin(I_r);
		VecAssemblyEnd(I_r);
		
		MatMult( mf, I_r, l );
		
		VecGetArray( l, &vals );
		
		nnz_local = nnz_local + M;
		if (tol) {
			/* Drop any small values */
			VecPointwiseMult( L, l,l );
			VecSqrt( L );
			VecGetArray( L, &_I_r );
			for( i=0; i<m; i++ ) {
				
#ifdef PETSC_USE_COMPLEX
				if( (PetscRealPart(_I_r[i])<tol) && (PetscImaginaryPart(_I_r[i])) ) {
					ridx[i] = -ridx[i];
					nnz_local--;
				}
#else
				if( _I_r[i] < tol ) {
					ridx[i] = -ridx[i];
					nnz_local--;
				}
#endif 
			}
			VecRestoreArray( L, &_I_r );
		}
		
		MatSetValues( A, m,ridx, 1,&c, vals, INSERT_VALUES );
		VecRestoreArray( l, &vals );
		
		if (tol) {
			/* undo negatives */
			for( i=0; i<m; i++ ) {
				if( ridx[i] < 0 ) {
					ridx[i] = -ridx[i];
				}
			}
		}
		
	}
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);	
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);	

	mierr=MPI_Allreduce( &nnz_local, &nnz, 1, MPI_INT, MPI_SUM, ((PetscObject)mf)->comm );CHKERRQ(mierr);
	PetscPrintf( ((PetscObject)mf)->comm, "MatFreeConstructExplicitOperator: Operator dimension %dx%d, nnz inserted %d \n", M,N,nnz );

	ierr=PetscFree( ridx );CHKERRQ(ierr);
	Stg_VecDestroy( &L );
	Stg_VecDestroy( &l );
	Stg_VecDestroy( &I_r );
	
	
	*_A = A;

	PetscFunctionReturn(0);
}

