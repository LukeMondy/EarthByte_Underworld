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
#include <assert.h>

#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>


#include "petscext.h"
#include "petscext_vec.h"
#include "petscext_mat.h"

/*
Mat GenerateMATAIJ( MPI_Comm comm, PetscInt M, PetscInt N, const char prefix[] )
{
	Mat A;
	
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetOptionsPrefix( A, prefix );
	MatSetSizes( A, PETSC_DECIDE,PETSC_DECIDE, M,N );
	MatSetType( A, MATAIJ );
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	return A;
}

void MatFillStride( Mat A, PetscScalar start_val, PetscScalar stride )
{
	PetscInt M,N, i,j;
	PetscScalar v;
	
	MatGetSize( A, &M, &N );
	v = start_val;
	for( i=0; i<M; i++ ) {
		for( j=0; j<N; j++ ) {
			MatSetValue( A, i,j, v, INSERT_VALUES );
			v = v + stride;
		}
	}
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}
*/


/*
[ A11 A12 xxx xxx ][phi]
[ A21 A22 A23 xxx ][P*]
[ A31 xxx  K   G  ][u]
[ xxx xxx  Gt  C  ][p]



*/


void test_magma_operator( PetscTruth view )
{
	Mat A;
	Mat Astokes, K,G,Gt,C;
	Mat Acomp, A11,A12,A21,A22;
	Mat A23, A31;
	PetscInt M,N,m,n;
	Vec x, phi,Pstar,u,p;
	Vec b, f_phi,f_Pstar,f,h;
	
	
	if( view ) PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	M = 10;
	N = 8;
	
	/* stokes objects */
	K = MATAIJNew( PETSC_COMM_WORLD, M,M, "K" );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G" );
	MatCreateSymTrans( PETSC_COMM_WORLD, G, &Gt );
	C = MATAIJNew( PETSC_COMM_WORLD, N,N, "C" );
	
	MatFillStride( K, 1.0, 2.0 );
	MatFillStride( G, 2.0, 4.4 );
	MatFillStride( C, 10.0, 0.0 );
	
	/* compaction objects */
	A11 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A11" );
	A12 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A12" );
	A21 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A21" );
	A22 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A22" );
	
	
	MatFillStride( A11, 1.0, 0.0 );
	MatFillStride( A12, 1.0, 1.0 );
	MatFillStride( A21, 1.0, 2.0 );
	MatFillStride( A22, 2.0, 0.0 );
	
	/* coupling terms */
	A23 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A23" );
	A31 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A31" );
	MatFillStride( A23, 10.0, 2.0 );
	MatFillStride( A31, 41.0, 55.0 );
	
	
	
	/* Assemble blocks */
	
	/* Global */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 4,4, 4,4 );
	MatSetType( A, "block" );
	
	MatSetOptionsPrefix( A, "A" );
	MatBlockSetValue( A, 0,0, A11,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatBlockSetValue( A, 2,2, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 2,3, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 3,2, Gt, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 3,3, C, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatBlockSetValue( A, 1,2, A23, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 2,0, A31, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	Stg_MatDestroy( K );		Stg_MatDestroy( G );
	Stg_MatDestroy( Gt );		Stg_MatDestroy( C );
	
	
	/* Get space for vectors */
	MatGetVecs( A, &x, &b );
	
	VecBlockGetSubVector( b, 0, &f_phi );
	VecBlockGetSubVector( b, 1, &f_Pstar );
	VecBlockGetSubVector( b, 2, &f );
	VecBlockGetSubVector( b, 3, &h );
	VecBlockRestoreSubVectors( b );
	
	VecBlockGetSubVector( x, 0, &phi );
	VecBlockGetSubVector( x, 1, &Pstar );
	VecBlockGetSubVector( x, 2, &u );
	VecBlockGetSubVector( x, 3, &p );
	VecBlockRestoreSubVectors( x );
	
	
	
	VecSet( f_phi, 		1.0 );
	VecSet( f_Pstar, 	2.0 );
	VecSet( f, 			3.0 );
	VecSet( h, 			4.0 );
	
	MatMult( A, b, x );
	if( view ) {
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
		VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	
	Stg_MatDestroy( A );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	
}



void test_magma_operator_2( PetscTruth view )
{
	Mat A;
	Mat Astokes, K,G,Gt,C;
	Mat Acomp, A11,A12,A21,A22;
	Mat A23, A31;
	PetscInt M,N,m,n;
	Vec x, phi,Pstar,u,p;
	Vec b, f_phi,f_Pstar,f,h;
	Vec x_comp, x_stokes;
	Vec b_comp, b_stokes;
	
	
	if( view ) PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	M = 10;
	N = 8;
	
	/* stokes objects */
	K = MATAIJNew( PETSC_COMM_WORLD, M,M, "K" );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G" );
	MatCreateSymTrans( PETSC_COMM_WORLD, G, &Gt );
	C = MATAIJNew( PETSC_COMM_WORLD, N,N, "C" );
	
	MatFillStride( K, 1.0, 2.0 );
	MatFillStride( G, 2.0, 4.4 );
	MatFillStride( C, 10.0, 0.0 );
	
	/* compaction objects */
	A11 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A11" );
	A12 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A12" );
	A21 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A21" );
//	A22 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A22" );
	
	
	MatFillStride( A11, 1.0, 0.0 );
	MatFillStride( A12, 1.0, 1.0 );
	MatFillStride( A21, 1.0, 2.0 );
//	MatFillStride( A22, 2.0, 0.0 );
	A22 = MATAIJIdentityNew( PETSC_COMM_WORLD, M );
	
	/* coupling terms */
	A23 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A23" );
	A31 = MATAIJNew( PETSC_COMM_WORLD, M,M, "A31" );
	MatFillStride( A23, 10.0, 2.0 );
	MatFillStride( A31, 41.0, 55.0 );
	
	
	
	/* Assemble blocks */
	
	/* Global */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 4,4, 4,4 );
	MatSetType( A, "block" );
	
	MatSetOptionsPrefix( A, "A" );
	MatBlockSetValue( A, 0,0, A11,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatBlockSetValue( A, 2,2, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 2,3, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 3,2, Gt, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 3,3, C, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatBlockSetValue( A, 1,2, A23, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A, 2,0, A31, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	/* 
	We hand off ownership of ALL sub matrices to A. Thus, A will release memory for sub matrices 
	*/
	Stg_MatDestroy( A11 );		Stg_MatDestroy( A12);
	Stg_MatDestroy( A21 );		Stg_MatDestroy( A22 );
	
	Stg_MatDestroy( K );		Stg_MatDestroy( G );
	Stg_MatDestroy( Gt );		Stg_MatDestroy( C );
	
	Stg_MatDestroy( A23 );		Stg_MatDestroy( A31 );
	
	/* Get space for vectors */
	MatGetVecs( A, &x, &b );
	
	/* 
	Memory for these vectors was created when we called MatGetVecs(), thus we will
	leave ownership with the global vectors x and b. We will NOT have to free the sub vectors
	*/
	VecBlockGetSubVector( b, 0, &f_phi );
	VecBlockGetSubVector( b, 1, &f_Pstar );
	VecBlockGetSubVector( b, 2, &f );
	VecBlockGetSubVector( b, 3, &h );
	VecBlockRestoreSubVectors( b );
	
	VecBlockGetSubVector( x, 0, &phi );
	VecBlockGetSubVector( x, 1, &Pstar );
	VecBlockGetSubVector( x, 2, &u );
	VecBlockGetSubVector( x, 3, &p );
	VecBlockRestoreSubVectors( x );
	
	
	/* Assemble compaction */
	MatCreate( PETSC_COMM_WORLD, &Acomp );
	MatSetSizes( Acomp, 2,2, 2,2 );
	MatSetType( Acomp, "block" );
	
	MatSetOptionsPrefix( Acomp, "Acomp" );
	MatBlockSetValue( Acomp, 0,0, A11,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Acomp, 0,1, A12,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Acomp, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( Acomp, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin( Acomp, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Acomp, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Assemble stokes */
	MatCreate( PETSC_COMM_WORLD, &Astokes );
	MatSetSizes( Astokes, 2,2, 2,2 );
	MatSetType( Astokes, "block" );
	
	MatSetOptionsPrefix( Astokes, "Astokes" );
	MatBlockSetValue( Astokes, 0,0, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Astokes, 0,1, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Astokes, 1,0, Gt, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( Astokes, 1,1, C, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin( Astokes, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Astokes, MAT_FINAL_ASSEMBLY);
	
	
	VecCreate( PETSC_COMM_WORLD, &x_comp );
	VecSetSizes( x_comp, 2,2 );
	VecSetType( x_comp, "block" );
	VecBlockSetValue( x_comp, 0, phi, INSERT_VALUES );
	VecBlockSetValue( x_comp, 1, Pstar, INSERT_VALUES );
	VecAssemblyBegin( x_comp );		VecAssemblyEnd( x_comp );
	
	VecDuplicate( x_comp, &b_comp );
	VecBlockSetValue( b_comp, 0, f_phi, INSERT_VALUES );
	VecBlockSetValue( b_comp, 1, f_Pstar, INSERT_VALUES );
	VecAssemblyBegin( b_comp );		VecAssemblyEnd( b_comp );
	
	
	VecCreate( PETSC_COMM_WORLD, &x_stokes );
	VecSetSizes( x_stokes, 2,2 );
	VecSetType( x_stokes, "block" );
	VecBlockSetValue( x_stokes, 0, u, INSERT_VALUES );
	VecBlockSetValue( x_stokes, 1, p, INSERT_VALUES );
	VecAssemblyBegin( x_stokes );		VecAssemblyEnd( x_stokes );
	
	VecDuplicate( x_stokes, &b_stokes );
	VecBlockSetValue( b_stokes, 0, f, INSERT_VALUES );
	VecBlockSetValue( b_stokes, 1, h, INSERT_VALUES );
	VecAssemblyBegin( b_stokes );		VecAssemblyEnd( b_stokes );
	
	
	
	/* Define value for vectors and do some multiplcations */
	VecSet( f_phi, 		1.0 );
	VecSet( f_Pstar, 	2.0 );
	VecSet( f, 			3.0 );
	VecSet( h, 			4.0 );
	
	MatMult( A, b, x );
	if( view ) {
		PetscPrintf( PETSC_COMM_WORLD, "[A]{b} = \n" );
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
		VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	MatMult( Astokes, b_stokes, x_stokes );
	if( view ) {
		PetscPrintf( PETSC_COMM_WORLD, "[Astokes]{b_stokes} = \n" );
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
		VecView( x_stokes, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	MatMult( Acomp, b_comp, x_comp );
	if( view ) {
		PetscPrintf( PETSC_COMM_WORLD, "[Acomp]{b_comp} = \n" );
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
		VecView( x_comp, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	if( view ) {
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
		MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	Stg_MatDestroy( Astokes );
	Stg_VecDestroy( & x_stokes );
	Stg_VecDestroy( & b_stokes );
	
	Stg_MatDestroy( Acomp );
	Stg_VecDestroy( & x_comp );
	Stg_VecDestroy( & b_comp );
	
	Stg_MatDestroy( A );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	
}






int main( int argc, char **args )
{
	int i;
	int BIG = 10000;
	PetscTruth test_mem, flg;
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	PetscExtInitialize();
	
	test_mem = PETSC_FALSE;
	PetscOptionsGetTruth( "", "-test_mem", &test_mem, &flg );
	
	if( test_mem ) {
		for( i=0; i<BIG; i++ ) {
		//	test_magma_operator(PETSC_FALSE);
			test_magma_operator_2(PETSC_FALSE);
		}
	}
	else {
		test_magma_operator(PETSC_TRUE);
		test_magma_operator_2(PETSC_TRUE);
	}
	
	
	
	PetscExtFinalize();
	PetscFinalize();
	return 0;
	
}
