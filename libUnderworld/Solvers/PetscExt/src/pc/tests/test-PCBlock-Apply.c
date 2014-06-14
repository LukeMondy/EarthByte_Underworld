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


#include "petscext_vec.h"
#include "petscext_mat.h"
#include "petscext_pc.h"


/*

Q = 
( 1/10   0      0   ) | (  0   3  )
(  0   1/100    0   ) | (  1   4  )
(  0     0   1/1000 ) | (  2   5  )
-----------------------------------
(   0    0      0   ) | ( 1/8  0  ) 
(   0    0      0   ) | (  0  1/9 )

X = ( 1 2 3 | 4 5 )^T
Y = Q^{-1} X 
Y = 
       -1340
      -21000
     -286000
          32
          45
*/
void test_pc_block_apply_U( void )
{
	int nu,np,i,j;
	Mat A11, A12, A22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	KSP ksp_Q;	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A11 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A11, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	/* Create + Setup A12 */
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A12, i,j, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A12 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A12, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	MatSetValue( A22, 0,0, (double)(1.0/8.0), INSERT_VALUES );
	MatSetValue( A22, 1,1, (double)(1.0/9.0), INSERT_VALUES );
	
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A22 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A22, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A12, &x2, &x1 );
	MatGetVecs( A12, &y2, &y1 );
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	for( i=0; i<nu; i++ ) {
		VecSetValue( x1, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( x1 );
	VecAssemblyEnd( x1 );
	for( i=0; i<np; i++ ) {
		VecSetValue( x2, i, (double)(1+i+nu), INSERT_VALUES );
	}
	VecAssemblyBegin( x2 );
	VecAssemblyEnd( x2 );

	
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( Y );
	VecAssemblyEnd( Y );
	
	
	
	
	
	/* Block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlockType( Q, PC_BLOCK_UPPER );
	
	/* configure */
	/* Q11 */
	PCBlockGetSubKSP( Q, 0, &ksp_Q );
	KSPGetPC( ksp_Q, &Q11 );
	PCSetType( Q11, "lu" );
	
	/* Q22 */
	PCBlockGetSubKSP( Q, 1, &ksp_Q );
	KSPGetPC( ksp_Q, &Q22 );
	PCSetType( Q22, "lu" );
	
	
	/* Apply preconditioner */
	PCApply( Q, X, Y );
	
	
//	VecView( Y, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockGetSubVector( Y, 0, &y1 );
	VecBlockGetSubVector( Y, 1, &y2 );
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A12 );	Stg_MatDestroy( A22 );
	Stg_MatDestroy( A );
	Stg_PCDestroy( Q );
	
}




/*



Q = 
( 1/10   0       0  ) | (  0   0 )
(  0   1/100     0  ) | (  0   0 )
(  0     0   1/1000 ) | (  0   0 )
----------------------------------
(  0     1       2   )| ( 1/8  0  ) 
(  3     4       5   )| (  0  1/9 )

X = ( 1 2 3 | 4 5 )^T
Y = Q^{-1} X 

Y =

          10
         200
        3000
      -49568
     -142425


*/


void test_pc_block_apply_L( void )
{
	int nu,np,i,j;
	Mat A11, A12, A21, A22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	KSP ksp_Q;	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	
	/* Create + Setup A12 */
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A12, i,j, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	MatTranspose( A12, MAT_INITIAL_MATRIX, &A21 );
	
	
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	MatSetValue( A22, 0,0, (double)(1.0/8.0), INSERT_VALUES );
	MatSetValue( A22, 1,1, (double)(1.0/9.0), INSERT_VALUES );
	
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A22 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A22, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A12, &x2, &x1 );
	MatGetVecs( A12, &y2, &y1 );
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	for( i=0; i<nu; i++ ) {
		VecSetValue( x1, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( x1 );
	VecAssemblyEnd( x1 );
	for( i=0; i<np; i++ ) {
		VecSetValue( x2, i, (double)(1+i+nu), INSERT_VALUES );
	}
	VecAssemblyBegin( x2 );
	VecAssemblyEnd( x2 );

	
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( Y );
	VecAssemblyEnd( Y );
	
	
	
	
	
	/* Block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlockType( Q, PC_BLOCK_LOWER );
	
	/* configure */
	/* Q11 */
	PCBlockGetSubKSP( Q, 0, &ksp_Q );
	KSPGetPC( ksp_Q, &Q11 );
	PCSetType( Q11, "lu" );
	
	/* Q22 */
	PCBlockGetSubKSP( Q, 1, &ksp_Q );
	KSPGetPC( ksp_Q, &Q22 );
	PCSetType( Q22, "lu" );
	
	
	/* Apply preconditioner */
	PCApply( Q, X, Y );
	
	
//	VecView( Y, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockGetSubVector( Y, 0, &y1 );
	VecBlockGetSubVector( Y, 1, &y2 );
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );	Stg_MatDestroy( A22 );
	Stg_MatDestroy( A );
	Stg_PCDestroy( Q );
	
}





#if 0

/*

A = 
( 1/10   0      0   ) | (  0   3  )
(  0   1/100    0   ) | (  1   4  )
(  0     0   1/1000 ) | (  2   5  )
-----------------------------------
(  0     1       2   )| (  0   0  ) 
(  3     4       5   )| (  0   0  )

Q11 = 
( 10   0      0   )
(  0   100    0   )
(  0     0   1000 ) 
-----------------------------------
S = ( 4100  10400 )
    ( 10400 26690 ) 


X = ( 1 2 3 | 4 5 )^T
Y = Q^{-1} X 

Y =
	-0.498818
	-2.3357
	 3.16785
	 0.623515
	 0.349961
Y =

  -0.49881796690307
  -2.33569739952719
   3.16784869976359
   0.62351457840820
   0.34996059889677



*/
void test_pc_block_apply_FULL( void )
{
	int nu,np,i,j;
	Mat A11, A12, A21, S22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	
	/* Create + Setup A12 */
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A12, i,j, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	
	/* Setup A21 */
	MatTranspose( A12, MAT_INITIAL_MATRIX, &A21 );
	
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, (nu+np),(nu+np), (nu+np),(nu+np) );
	MatSetType( A, "block" );
	MatBlock_SetOperators( A,
			A11, A12, A21, PETSC_NULL, 
			PETSC_NULL, PETSC_NULL );
	
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A12, &x2, &x1 );
	MatGetVecs( A12, &y2, &y1 );
	MatBlock_CreateGlobalVector( A, &X );
	VecDuplicate( X, &Y );
	
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	VecZeroEntries( X );
	for( i=0; i<nu+np; i++ ) {
		VecSetValue( X, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	
	
	/* Create pc, Q11 */
	PCCreate( PETSC_COMM_WORLD, &Q11 );
	PCSetOperators( Q11, A11, A11, SAME_NONZERO_PATTERN );
	PCSetType( Q11, "lu" );
	
	
	/* Create pc, Q22 */
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &S22 );
	MatSetSizes( S22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( S22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( S22, np, PETSC_NULL );
	
	MatSetValue( S22, 0,0, (double)(4100), INSERT_VALUES );
	MatSetValue( S22, 0,1, (double)(10400), INSERT_VALUES );
	MatSetValue( S22, 1,0, (double)(10400), INSERT_VALUES );
	MatSetValue( S22, 1,1, (double)(26690), INSERT_VALUES );
	


	MatAssemblyBegin( S22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( S22, MAT_FINAL_ASSEMBLY);
	
	
	PCCreate( PETSC_COMM_WORLD, &Q22 );
	PCSetOperators( Q22, S22, S22, SAME_NONZERO_PATTERN );
	PCSetType( Q22, "lu" );
	
	
	/* Setup the block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlocks( Q, Q11, A12, A21, Q22 );
	PCBlockSetBlockType( Q, PC_BLOCK_FULL );
	
	
	/* Apply preconditioner */
	PCApply( Q, X, Y );
	
	VecView( Y, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_VecDestroy( & x1 );	Stg_VecDestroy( & x2 );
	
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );	Stg_MatDestroy( S22 );
	Stg_MatDestroy( A );
	
	
	Stg_PCDestroy( Q11 );	Stg_PCDestroy( Q22 );
	Stg_PCDestroy( Q );
	
}

#endif


/*
Should be same as test_pc_block_apply_U()
*/
void test_pc_block_apply_LTrans( void )
{
	int nu,np,i,j;
	Mat A11, A12, A21, A22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	KSP ksp_Q;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	
	/* Create + Setup A12 */
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A12, i,j, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	MatTranspose( A12, MAT_INITIAL_MATRIX, &A21 );
	
	
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	MatSetValue( A22, 0,0, (double)(1.0/8.0), INSERT_VALUES );
	MatSetValue( A22, 1,1, (double)(1.0/9.0), INSERT_VALUES );
	
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A22 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A22, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A12, &x2, &x1 );
	MatGetVecs( A12, &y2, &y1 );
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	for( i=0; i<nu; i++ ) {
		VecSetValue( x1, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( x1 );
	VecAssemblyEnd( x1 );
	for( i=0; i<np; i++ ) {
		VecSetValue( x2, i, (double)(1+i+nu), INSERT_VALUES );
	}
	VecAssemblyBegin( x2 );
	VecAssemblyEnd( x2 );

	
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( Y );
	VecAssemblyEnd( Y );
	
	
	
	
	
	/* Block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlockType( Q, PC_BLOCK_LOWER );
	
	/* configure */
	/* Q11 */
	PCBlockGetSubKSP( Q, 0, &ksp_Q );
	KSPGetPC( ksp_Q, &Q11 );
	PCSetType( Q11, "lu" );
	
	/* Q22 */
	PCBlockGetSubKSP( Q, 1, &ksp_Q );
	KSPGetPC( ksp_Q, &Q22 );
	PCSetType( Q22, "lu" );
	
	
	/* Apply preconditioner */
	PCApplyTranspose( Q, X, Y );
	
	
//	VecView( Y, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockGetSubVector( Y, 0, &y1 );
	VecBlockGetSubVector( Y, 1, &y2 );
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );	Stg_MatDestroy( A22 );
	Stg_MatDestroy( A );
	Stg_PCDestroy( Q );
	
}



/*
Should be same as test_pc_block_apply_L()
*/
void test_pc_block_apply_UTrans( void )
{
	int nu,np,i,j;
	Mat A11, A12, A22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	KSP ksp_Q;	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A11 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A11, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	/* Create + Setup A12 */
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A12, i,j, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A12 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A12, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	MatSetValue( A22, 0,0, (double)(1.0/8.0), INSERT_VALUES );
	MatSetValue( A22, 1,1, (double)(1.0/9.0), INSERT_VALUES );
	
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A22 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A22, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A12, &x2, &x1 );
	MatGetVecs( A12, &y2, &y1 );
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	for( i=0; i<nu; i++ ) {
		VecSetValue( x1, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( x1 );
	VecAssemblyEnd( x1 );
	for( i=0; i<np; i++ ) {
		VecSetValue( x2, i, (double)(1+i+nu), INSERT_VALUES );
	}
	VecAssemblyBegin( x2 );
	VecAssemblyEnd( x2 );

	
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( Y );
	VecAssemblyEnd( Y );
	
	
	
	
	
	/* Block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlockType( Q, PC_BLOCK_UPPER );
	
	/* configure */
	/* Q11 */
	PCBlockGetSubKSP( Q, 0, &ksp_Q );
	KSPGetPC( ksp_Q, &Q11 );
	PCSetType( Q11, "lu" );
	
	/* Q22 */
	PCBlockGetSubKSP( Q, 1, &ksp_Q );
	KSPGetPC( ksp_Q, &Q22 );
	PCSetType( Q22, "lu" );
	
	
	/* Apply preconditioner */
	PCApplyTranspose( Q, X, Y );
	
	
//	VecView( Y, PETSC_VIEWER_STDOUT_WORLD );
	VecBlockGetSubVector( Y, 0, &y1 );
	VecBlockGetSubVector( Y, 1, &y2 );
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A12 );	Stg_MatDestroy( A22 );
	Stg_MatDestroy( A );
	Stg_PCDestroy( Q );
	
}


/*


A21 = ( 0 1 2 )
      ( 3 4 5 )


Q = 
( 1/10   0      0   ) | (         )
(  0   1/100    0   ) | (  A21^T  )
(  0     0   1/1000 ) | (         )
-----------------------------------
(   0    0      0   ) | ( 1/8  0  ) 
(   0    0      0   ) | (  0  1/9 )

X = ( 1 2 3 | 4 5 )^T
Y = Q^{-1} X 
Y = 
       -1340
      -21000
     -286000
          32
          45
*/


void test_pc_block_apply_USymTrans( void )
{
	int nu,np,i,j;
	Mat A11, A21Trans, A21, A22, A;
	Vec x1,x2, X,Y, y1,y2;
	PC Q11, Q22, Q;
	KSP ksp_Q;	
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	nu = 3;
	np = 2;
	
	
	/* Create + Setup A11 */
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, nu, nu );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, nu, PETSC_NULL );
	
	MatSetValue( A11, 0,0, (double)(1.0/10.0), INSERT_VALUES );
	MatSetValue( A11, 1,1, (double)(1.0/100.0), INSERT_VALUES );
	MatSetValue( A11, 2,2, (double)(1.0/1000.0), INSERT_VALUES );
	
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	
	/* Create + Setup A21 */
	MatCreate( PETSC_COMM_WORLD, &A21 );
	MatSetSizes( A21, PETSC_DECIDE, PETSC_DECIDE, np, nu );
	MatSetType( A21, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A21, nu, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( A21, j,i, (double)(i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A21, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A21, MAT_FINAL_ASSEMBLY);
	
	/* Create symbolic A12 */
	MatCreateSymTrans( PETSC_COMM_WORLD, A21, &A21Trans );
	
	
	/* create a matrix for the (2,2) block */
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, np, np );
	MatSetType( A22, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A22, np, PETSC_NULL );
	
	MatSetValue( A22, 0,0, (double)(1.0/8.0), INSERT_VALUES );
	MatSetValue( A22, 1,1, (double)(1.0/9.0), INSERT_VALUES );
	
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	/*
	PetscPrintf( PETSC_COMM_WORLD, "A22 = \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A22, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Create + Setup A */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A21Trans, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create + Setup vectors */
	MatGetVecs( A21Trans, &x2, &x1 );
	MatGetVecs( A21Trans, &y2, &y1 );
	
	/* Setup vectors X */
	VecZeroEntries( x1 );
	VecZeroEntries( x2 );
	for( i=0; i<nu; i++ ) {
		VecSetValue( x1, i, (double)(1+i), INSERT_VALUES );
	}
	VecAssemblyBegin( x1 );
	VecAssemblyEnd( x1 );
	for( i=0; i<np; i++ ) {
		VecSetValue( x2, i, (double)(1+i+nu), INSERT_VALUES );
	}
	VecAssemblyBegin( x2 );
	VecAssemblyEnd( x2 );

	
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2,2 );
	VecSetType( X, "block" );
	
	VecBlockSetValue( X, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( X, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( X );
	VecAssemblyEnd( X );
	
	VecCreate( PETSC_COMM_WORLD, &Y );
	VecSetSizes( Y, 2,2 );
	VecSetType( Y, "block" );
	
	VecBlockSetValue( Y, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( Y, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( Y );
	VecAssemblyEnd( Y );
	
	
	
	
	
	/* Block pc */
	PCCreate( PETSC_COMM_WORLD, &Q );
	PCSetOperators( Q, A, A, SAME_NONZERO_PATTERN );
	PCSetType( Q, "block" );
	PCBlockSetBlockType( Q, PC_BLOCK_UPPER );
	
	/* configure */
	/* Q11 */
	PCBlockGetSubKSP( Q, 0, &ksp_Q );
	KSPGetPC( ksp_Q, &Q11 );
	PCSetType( Q11, "jacobi" );
	
	/* Q22 */
	PCBlockGetSubKSP( Q, 1, &ksp_Q );
	KSPGetPC( ksp_Q, &Q22 );
	PCSetType( Q22, "jacobi" );
	
	
	/* Apply preconditioner */
	PCApply( Q, X, Y );
	
	
	VecBlockGetSubVector( Y, 0, &y1 );
	VecBlockGetSubVector( Y, 1, &y2 );
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	Stg_VecDestroy( & X );	Stg_VecDestroy( & Y );
	Stg_MatDestroy( A11 );	Stg_MatDestroy( A21Trans );
	Stg_MatDestroy( A21 );	Stg_MatDestroy( A22 );
	Stg_MatDestroy( A );
	Stg_PCDestroy( Q );
	
}




int main( int argc, char **args )
{
	int i;
	int BIG = 1;
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	PetscExtPCRegisterAll();
	
	
	for( i=0; i<BIG; i++ ) {
		test_pc_block_apply_U();
		test_pc_block_apply_L();
		
		test_pc_block_apply_LTrans();
		test_pc_block_apply_UTrans();
		
		test_pc_block_apply_USymTrans();
		
		/* Not mathematically formulated correctly yet */
//		test_pc_block_apply_FULL();
		
	}
	
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	PetscExtPCRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
	
}
