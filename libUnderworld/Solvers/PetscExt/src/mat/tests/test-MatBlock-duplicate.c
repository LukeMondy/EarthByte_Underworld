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


void test_MatBlock_duplicate( void )
{
	Mat A11,A12,A21,A22, A,B;
	Vec diag;
	int n, np, i,j;
	const MatType type;
	
	
//	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
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
	
	// A22
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	for( i=0; i<np; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A22, i,j, (double)(i+j*n), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	
	
	/* Create block matrix */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	MatGetType( A, &type );
	
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	
	
	// A12
	MatBlockCreateSubMatrix( A, MATAIJ, 0,1, &A12 );
	// A21
	MatBlockCreateSubMatrix( A, MATAIJ, 1,0, &A21 );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A12, i,j, 1.0, INSERT_VALUES );
			MatSetValue( A21, j,i, 1.0, INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin( A21, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A21, MAT_FINAL_ASSEMBLY);
	
	
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	Stg_MatDestroy( A11 );		Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );		Stg_MatDestroy( A22 );
	
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
//	printf("A \n");
//	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	MatDuplicate( A, MAT_COPY_VALUES, &B );
//	printf("B \n");
//	MatView( B, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_MatDestroy( A );
	Stg_MatDestroy( B );
}



int main( int argc, char **args )
{
	int i;
	int BIG=1;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	for( i=0; i<10000; i++ ) {
		test_MatBlock_duplicate();
	}
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
