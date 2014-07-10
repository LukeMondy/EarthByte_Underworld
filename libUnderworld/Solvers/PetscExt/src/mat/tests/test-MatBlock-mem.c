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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>

#include "petscext_vec.h"
#include "petscext_mat.h"


void test_solve( void )
{
	Mat A11, A12,A21, A;
	KSP ksp;
	PC pc;
	Vec b,x , f,h, diag, x1,x2;
	int n, np, i,j;
	const MatType type;
	
	
	
	n = 3;
	np = 2;
	
	/* Create matrices */
	// A11
	VecCreate( PETSC_COMM_WORLD, &diag );
	VecSetSizes( diag, PETSC_DECIDE, n );
	VecSetType( diag, VECSEQ );
	
	VecSet( diag, (1.0/10.0) ); /* so inverse = diag(10) */
	
	/* As a test, create a diagonal matrix for A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, n, n );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, n, PETSC_NULL );
	MatDiagonalSet( A11, diag, INSERT_VALUES );
	
	Stg_VecDestroy( & diag );
	
	// A12
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, n, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A12, i,j, (double)(i+j*n), INSERT_VALUES );
		}
	}
	MatSetValue( A12, 2,1, (double)(4), INSERT_VALUES );
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	// A21
	MatTranspose( A12, MAT_INITIAL_MATRIX, &A21 );
	
	
	
	/* Create block matrix */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes_Block( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	MatGetType( A, &type );
	
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	
	/* Create vectors */
	MatGetVecs( A12, &h, &f );
	
	VecSet( f, 1.0 );
	VecSet( h, 0.0 );
	
	
	/* Create block vector */
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, 2, 2 );
	VecSetType( b, "block" );
	
	VecBlockSetValue( b, 0, f, INSERT_VALUES );
	VecBlockSetValue( b, 1, h, INSERT_VALUES );
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	
	VecDuplicate( b, &x );
	
	
	KSPCreate( PETSC_COMM_WORLD, &ksp );
	Stg_KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPSetType( ksp, "gmres" );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, "none" );
	KSPSetFromOptions( ksp );
	
	KSPSolve( ksp, b, x );
	
	
	VecBlockGetSubVector( x, 0, &x1 );
	VecBlockGetSubVector( x, 1, &x2 );
	VecBlockRestoreSubVectors( x );
	
	KSPDestroy( & ksp );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	Stg_MatDestroy( A11 );
	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );
	Stg_VecDestroy( & f );
	Stg_VecDestroy( & h );
	
	Stg_MatDestroy( A );
	
}



int main( int argc, char **args )
{
	int i;
	int BIG=10000;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	for( i=0; i<BIG; i++ ) {
		test_solve();
	}
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
