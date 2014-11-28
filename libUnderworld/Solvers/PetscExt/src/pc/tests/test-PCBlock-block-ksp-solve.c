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


/*
Suppose we wish to solve the block system,
  [A] = diag( [a11], a[22], [a33] )
We could 
 i)  Apply a Krylov method to the entire system
 ii) As the blocks are completely decoupled, we could apply a Krylov method to each block

I don't yet know if there is a serious disadvatage to using choice (i). Hence I wish to ensure
that the block framework can deliver results using the second options.

In this test, we experiment with the second option.

*/



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
#include "petscext_pc.h"




void test_block_solve( void )
{
	Mat A11, A22, Amat;
	PC pc, pc_diag;
	KSP ksp_diag;
	KSP ksp;
	Vec x1,x2,y1,y2, bX, bY;
	int n, n22, i,j;
	PetscInt nb;
	
	
	
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
		MatSetValue( A22, i,i, (double)(450.0 + i), INSERT_VALUES );
	}
	MatAssemblyBegin( A22, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A22, MAT_FINAL_ASSEMBLY);
	
	
	
	/* Create vectors */
	MatGetVecs( A11, PETSC_NULL, &x1 );
	MatGetVecs( A22, PETSC_NULL, &x2 );
	
	VecDuplicate( x1, &y1 );
	VecDuplicate( x2, &y2 );
	
	
	VecSet( x1, 10.0 );
	VecSet( x2, 55.0 );
	
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
	
	
	
	/* Solve A y = x */
	KSPCreate( PETSC_COMM_WORLD, &ksp );
	Stg_KSPSetOperators( ksp, Amat, Amat, SAME_NONZERO_PATTERN );
	KSPSetType( ksp, "preonly" );
	
	KSPGetPC( ksp, &pc );
	PCSetType( pc, "block" );
	
	MatGetSize( Amat, &nb, &nb );
	for( i=0; i<nb; i++ ) {
		PCBlockGetSubKSP( pc, i, &ksp_diag );
		KSPSetType( ksp_diag, "cg" );
	}
	
	
	
	KSPSolve( ksp, bX, bY );
	
	/*
	PetscPrintf( PETSC_COMM_WORLD, "y1 = \n");
	VecView( y1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "y2 = \n");
	VecView( y2, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	KSPDestroy( & ksp );
	
	
	Stg_VecDestroy( & bX );
	Stg_VecDestroy( & bY );
	Stg_MatDestroy( Amat );
	
	Stg_MatDestroy( A11 );
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
	
	
	test_block_solve();
	
	PetscOptionsGetTruth( PETSC_NULL, "-leak", &val, &flg );
	if( flg == PETSC_TRUE ) {
		PetscPrintf( PETSC_COMM_WORLD, "  Performing leak test \n");
		for( i=0; i<BIG; i++ ) {
			
			test_block_solve();
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
