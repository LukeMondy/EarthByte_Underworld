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


void test_MatBlock_merge( PetscInt FAC, PetscTruth view_merged_mat )
{
	Mat A11,A12,A22, A,B;
	Vec diag;
	int n, np, i,j;
	const MatType type;
	

	
	n = 3*FAC;
	np = 2*FAC;
	
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
			MatSetValue( A22, i,j, (double)(j+i*np + 1.0), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	
	
	/* Create block matrix */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes_Block( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	MatGetType( A, &type );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	printf("Forming A12 \n");
	// A12
	MatBlockCreateSubMatrix( A, MATAIJ, 0,1, &A12 );
        MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A12, i,j, 100.0*(j+i*np+1.0), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	MatBlockSetValue( A, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	if( view_merged_mat == PETSC_TRUE ) {
          PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
          printf("A11 \n");
          MatView( A11, PETSC_VIEWER_STDOUT_WORLD );
          printf("A12 \n");
          MatView( A12, PETSC_VIEWER_STDOUT_WORLD );
          printf("A22 \n");
          MatView( A22, PETSC_VIEWER_STDOUT_WORLD );	
	}

	Stg_MatDestroy( A11 );		Stg_MatDestroy( A12 );
	Stg_MatDestroy( A22 );
	

	printf("Merging B <INSERT>\n");
	/*
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	printf("A \n");
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	B = PETSC_NULL;
	MatBlockMergeSubBlocks( A, INSERT_VALUES, MATSEQAIJ, &B );

	if(view_merged_mat==PETSC_TRUE) {
	  printf("B \n");
	  PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	  MatView( B, PETSC_VIEWER_STDOUT_WORLD );
	}	


	printf("Merging B <ADD>\n");
        MatBlockMergeSubBlocks( A, ADD_VALUES, MATAIJ, &B );

	if(view_merged_mat==PETSC_TRUE) {
          printf("B <= B + B\n");
          PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
          MatView( B, PETSC_VIEWER_STDOUT_WORLD );
	}

	{
	  Mat AA;
	
	  printf("MatConvert BLOCK => MATAIJ\n");
	  MatConvert( A, MATAIJ, MAT_INITIAL_MATRIX, &AA );

	/*
	// CANNOT DO THIS. POINTER FOR A MUST MATCH AA //
          MatConvert( A, MATAIJ, MAT_REUSE_MATRIX, &AA );
          printf("AA <reuse>\n");
          PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
          MatView( AA, PETSC_VIEWER_STDOUT_WORLD );
	*/


	  if(view_merged_mat==PETSC_TRUE) {
            printf("AA \n");
            PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
            MatView( AA, PETSC_VIEWER_STDOUT_WORLD );       
	  }

	  Stg_MatDestroy(AA);
	}


	Stg_MatDestroy( A );
	Stg_MatDestroy( B );
}



int main( int argc, char **args )
{
	int i;
	int BIG=1;
	PetscInt FAC;
	PetscTruth view,flg;	

	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	view = PETSC_FALSE;
	PetscOptionsGetTruth( PETSC_NULL, "-view_merged_mat", &view, &flg );

	FAC = 10;
	PetscOptionsGetInt( PETSC_NULL, "-FAC", &FAC, &flg );

	for( i=0; i<BIG; i++ ) {
		test_MatBlock_merge(FAC,view);
	}
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
