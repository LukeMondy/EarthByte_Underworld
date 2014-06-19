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


#include "petscext.h"
#include "petscext_vec.h"
#include "petscext_mat.h"
#include "petscext_pc.h"


void test_singleton( void )
{
	Mat A11, Amat;
	PC pc;
	Vec x,y, bX, bY;
	int n, i,j;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	n = 10;
	
	/* Create matrices */
	// A11
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, n, n );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, n, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		/*
		for( j=0; j<n; j++ ) {
			MatSetValue( A11, i,j, (double)(10.0 + i), INSERT_VALUES );
		}
		*/
		MatSetValue( A11, i,i, (double)(10.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	/* Create vectors */
	MatGetVecs( A11, &x, &y );
	
	VecSet( x, 1.0 );
	VecSet( y, 0.0 );
	
	
	/* Create a block matrix */
	MatCreate( PETSC_COMM_WORLD, &Amat );
	MatSetSizes_Block( Amat, 1,1, 1,1 );
	MatSetType( Amat, "block" );
	MatBlockSetValue( Amat, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY);
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &bX );
	VecSetSizes( bX, 1,1 );
	VecSetType( bX, "block" );
	VecBlockSetValue( bX, 0, x, INSERT_VALUES );	Stg_VecDestroy( & x );
	VecAssemblyBegin( bX );
	VecAssemblyEnd( bX );
	
	VecCreate( PETSC_COMM_WORLD, &bY );
	VecSetSizes( bY, 1,1 );
	VecSetType( bY, "block" );
	VecBlockSetValue( bY, 0, y, INSERT_VALUES );	Stg_VecDestroy( & y );
	VecAssemblyBegin( bY );
	VecAssemblyEnd( bY );
	
	
	PCCreate( PETSC_COMM_WORLD, &pc );
	PCSetOperators( pc, Amat, Amat, SAME_NONZERO_PATTERN );
	PCSetType( pc, "block" );
	
	
	PCApply( pc, bX, bY );
	//PCView( pc, PETSC_VIEWER_STDOUT_WORLD );
	
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_PCDestroy( pc );
	
	Stg_VecDestroy( & bX );
	Stg_VecDestroy( & bY );
	
	Stg_MatDestroy( Amat );
	Stg_MatDestroy( A11 );
}


void test_2block_DIAG( void )
{
	Mat A11, A22, Amat;
	PC pc, _pc_i;
	Vec x1,x2,y1,y2, bX, bY;
	int n, n22, i,j;
	const char *bt;
	const PCType pc_type;
	KSP _ksp_i;	

	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	n = 10;
	n22 = 12;
	
	/* Create matrices */
	// A11
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, n, n );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, n, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		MatSetValue( A11, i,i, (double)(10.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	// A22
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, n22, n22 );
	MatSetType( A22, MATSEQAIJ );
	
	for( i=0; i<n22; i++ ) {
		MatSetValue( A22, i,i, (double)(10.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create vectors */
	MatGetVecs( A11, PETSC_NULL, &x1 );
	MatGetVecs( A22, PETSC_NULL, &x2 );
	
	VecDuplicate( x1, &y1 );
	VecDuplicate( x2, &y2 );
	
	
	VecSet( x1, 1.0 );
	VecSet( x2, 1.0 );
	
	VecSet( y1, 0.0 );
	VecSet( y2, 0.0 );
	
	
	/* Create a block matrix */
	MatCreate( PETSC_COMM_WORLD, &Amat );
	MatSetSizes_Block( Amat, 2,2, 2,2 );
	MatSetType( Amat, "block" );
	MatBlockSetValue( Amat, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Amat, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY);
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &bX );
	VecSetSizes( bX, 2,2 );
	VecSetType( bX, "block" );
	
	VecBlockSetValue( bX, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( bX, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( bX );
	VecAssemblyEnd( bX );
	
	VecCreate( PETSC_COMM_WORLD, &bY );
	VecSetSizes( bY, 2,2 );
	VecSetType( bY, "block" );
	
	VecBlockSetValue( bY, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( bY, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( bY );
	VecAssemblyEnd( bY );
	
	
	PCCreate( PETSC_COMM_WORLD, &pc );
	PCSetOperators( pc, Amat, Amat, SAME_NONZERO_PATTERN );
	PCSetType( pc, "block" );
	
	PCBlockGetSubKSP( pc, 0, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCSetType( _pc_i, "lu" );
	PCBlockGetSubKSP( pc, 1, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCSetType( _pc_i, "jacobi" );
	
	
	PCApply( pc, bX, bY );
	//PCView( pc, PETSC_VIEWER_STDOUT_WORLD );
	
	
	PCBlockGetBlockType( pc, &bt );
	PetscPrintf( PETSC_COMM_WORLD, "block_type: %s \n", bt );
	PCBlockGetSubKSP( pc, 0, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCGetType( _pc_i, &pc_type );
	PetscPrintf( PETSC_COMM_WORLD, "(0,0) pc_type: %s \n", pc_type );
	
	PCBlockGetSubKSP( pc, 1, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCGetType( _pc_i, &pc_type );
	PetscPrintf( PETSC_COMM_WORLD, "(1,1) pc_type: %s \n", pc_type );
	
	
	PetscPrintf( PETSC_COMM_WORLD, "y1 = \n");
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "y2 = \n");
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_PCDestroy( pc );
	
	Stg_VecDestroy( & bX );
	Stg_VecDestroy( & bY );
	Stg_MatDestroy( Amat );
	
	
	Stg_MatDestroy( A11 );
	Stg_MatDestroy( A22 );
}



void test_2block_UPPER( void )
{
	Mat A11, A22, A12, Amat;
	PC pc, _pc_i;
	Vec x1,x2,y1,y2, bX, bY;
	int n, n22, i,j;
	const char *bt;
	const PCType pc_type;
	KSP _ksp_i;	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	
	n = 10;
	n22 = 12;
	
	/* Create matrices */
	// A11
	MatCreate( PETSC_COMM_WORLD, &A11 );
	MatSetSizes( A11, PETSC_DECIDE, PETSC_DECIDE, n, n );
	MatSetType( A11, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A11, n, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		MatSetValue( A11, i,i, (double)(10.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A11, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A11, MAT_FINAL_ASSEMBLY);
	
	// A22
	MatCreate( PETSC_COMM_WORLD, &A22 );
	MatSetSizes( A22, PETSC_DECIDE, PETSC_DECIDE, n22, n22 );
	MatSetType( A22, MATSEQAIJ );
	
	for( i=0; i<n22; i++ ) {
		MatSetValue( A22, i,i, (double)(10.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	
	
	// A12
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, n, n22 );
	MatSetType( A12, MATSEQAIJ );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<n22; j++ ) {
			MatSetValue( A12, i,j, (double)(i+j*n), INSERT_VALUES );
	}}
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	
	
	
	/* Create vectors */
	MatGetVecs( A11, PETSC_NULL, &x1 );
	MatGetVecs( A22, PETSC_NULL, &x2 );
	
	VecDuplicate( x1, &y1 );
	VecDuplicate( x2, &y2 );
	
	
	VecSet( x1, 1.0 );
	VecSet( x2, 1.0 );
	
	VecSet( y1, 0.0 );
	VecSet( y2, 0.0 );
	
	
	/* Create a block matrix */
	MatCreate( PETSC_COMM_WORLD, &Amat );
	MatSetSizes_Block( Amat, 2,2, 2,2 );
	MatSetType( Amat, "block" );
	
	MatBlockSetValue( Amat, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Amat, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( Amat, 1,1, A22, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	
	MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY);
	
	/* Create block vectors */
	VecCreate( PETSC_COMM_WORLD, &bX );
	VecSetSizes( bX, 2,2 );
	VecSetType( bX, "block" );
	
	VecBlockSetValue( bX, 0, x1, INSERT_VALUES );	Stg_VecDestroy( & x1 );
	VecBlockSetValue( bX, 1, x2, INSERT_VALUES );	Stg_VecDestroy( & x2 );
	
	VecAssemblyBegin( bX );
	VecAssemblyEnd( bX );
	
	VecCreate( PETSC_COMM_WORLD, &bY );
	VecSetSizes( bY, 2,2 );
	VecSetType( bY, "block" );
	
	VecBlockSetValue( bY, 0, y1, INSERT_VALUES );	Stg_VecDestroy( & y1 );
	VecBlockSetValue( bY, 1, y2, INSERT_VALUES );	Stg_VecDestroy( & y2 );
	
	VecAssemblyBegin( bY );
	VecAssemblyEnd( bY );
	
	
	PCCreate( PETSC_COMM_WORLD, &pc );
	PCSetOperators( pc, Amat, Amat, SAME_NONZERO_PATTERN );
	PCSetType( pc, "block" );
	PCBlockSetBlockType( pc, PC_BLOCK_UPPER );
	
	PCApply( pc, bX, bY );
//	PCView( pc, PETSC_VIEWER_STDOUT_WORLD );
	
	PCBlockGetBlockType( pc, &bt );
	PetscPrintf( PETSC_COMM_WORLD, "block_type: %s \n", bt );
	PCBlockGetSubKSP( pc, 0, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCGetType( _pc_i, &pc_type );
	PetscPrintf( PETSC_COMM_WORLD, "(0,0) pc_type: %s \n", pc_type );
	
	PCBlockGetSubKSP( pc, 1, &_ksp_i );
	KSPGetPC( _ksp_i, &_pc_i );
	PCGetType( _pc_i, &pc_type );
	PetscPrintf( PETSC_COMM_WORLD, "(1,1) pc_type: %s \n", pc_type );
	
	PetscPrintf( PETSC_COMM_WORLD, "y1 = \n");
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "y2 = \n");
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );

	
	Stg_PCDestroy( pc );
	
	
	Stg_VecDestroy( & bX );
	Stg_VecDestroy( & bY );
	Stg_MatDestroy( Amat );
	
	
	Stg_MatDestroy( A11 );
	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A22 );
}



int main( int argc, char **args )
{
	int i;
	int BIG=10000;
	PetscTruth flg, val;
	
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	PetscExtPCRegisterAll();
	
	
	
	test_singleton();
	test_2block_DIAG();
	test_2block_UPPER();
//	test_2block_LOWER();
	
	PetscOptionsGetTruth( PETSC_NULL, "-leak", &val, &flg );
	if( flg == PETSC_TRUE ) {
		PetscPrintf( PETSC_COMM_WORLD, "  Performing leak test \n");
		for( i=0; i<BIG; i++ ) {
			
		//	test_singleton();
			test_2block_DIAG();
			
			if( i % (BIG/10) == 0 ) {
				PetscPrintf( PETSC_COMM_WORLD, "     Tests conducted: %d/%d \n", i,BIG );
			}
		}
	}
	
	
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	PetscExtPCRegisterDestroyAll();
	PetscFinalize();
	return 0;
}
