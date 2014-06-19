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


/*

Q = 
(  1   4  )
(  2   5  )
(  3   6  )

*/
void test_sym_trans_Apply( void )
{
	int nu,np,i,j;
	Mat Q, symQt, Qt;
	Vec right, left;
	
	
	nu = 100;
	np = 70;
	
	
	/* Create + Setup Q */
	MatCreate( PETSC_COMM_WORLD, &Q );
	MatSetSizes( Q, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( Q, MATAIJ );
	MatSeqAIJSetPreallocation( Q, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( Q, i,j, (double)(1.0 + i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( Q, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Q, MAT_FINAL_ASSEMBLY);
	
	MatTranspose( Q, MAT_INITIAL_MATRIX, &Qt );
	
	
	MatGetVecs( Qt, &right, &left );
	VecSet( right, 1.0 );
	
	MatMult( Qt, right, left );
	
	MatCreateSymTrans( PETSC_COMM_WORLD, Q, &symQt );
	MatMult( symQt, right, left );
	
	
	Stg_VecDestroy( & right );	Stg_VecDestroy( & left );
	Stg_MatDestroy( Q );
	Stg_MatDestroy( Qt );
	Stg_MatDestroy( symQt );
}


void test_sym_trans_ApplyTrans( void )
{
	int nu,np,i,j;
	Mat Q, symQt, Qt;
	Vec right, left;
	
	
	nu = 100;
	np = 70;
	
	
	/* Create + Setup Q */
	MatCreate( PETSC_COMM_WORLD, &Q );
	MatSetSizes( Q, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( Q, MATAIJ );
	MatSeqAIJSetPreallocation( Q, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( Q, i,j, (double)(1.0 + i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( Q, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Q, MAT_FINAL_ASSEMBLY);
	
	MatGetVecs( Q, &right, &left );
	VecSet( right, 1.0 );
	
	MatMult( Q, right, left );
	
	MatCreateSymTrans( PETSC_COMM_WORLD, Q, &symQt );
	MatMultTranspose( symQt, right, left );
	
	Stg_VecDestroy( & right );	Stg_VecDestroy( & left );
	Stg_MatDestroy( Q );
	Stg_MatDestroy( symQt );
}



void test_sym_trans_GetExplicitOp( void )
{
	int nu,np,i,j;
	Mat Q, symQt;
	Mat exp_op;
	
	nu = 100;
	np = 70;
	
	
	/* Create + Setup Q */
	MatCreate( PETSC_COMM_WORLD, &Q );
	MatSetSizes( Q, PETSC_DECIDE, PETSC_DECIDE, nu, np );
	MatSetType( Q, MATAIJ );
	MatSeqAIJSetPreallocation( Q, np, PETSC_NULL );
	
	for( j=0; j<np; j++ ) {
		for( i=0; i<nu; i++ ) {
			MatSetValue( Q, i,j, (double)(1.0 + i+j*nu), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( Q, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( Q, MAT_FINAL_ASSEMBLY);
	
	
	
	
	MatCreateSymTrans( PETSC_COMM_WORLD, Q, &symQt );
	MatSymTransGetExplicitOperator( symQt, &exp_op );
	
	Stg_MatDestroy( Q );
	Stg_MatDestroy( symQt );
	Stg_MatDestroy( exp_op );
}



int main( int argc, char **args )
{
	int i;
	int BIG = 1000;
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	
	for( i=0; i<BIG; i++ ) {
		test_sym_trans_Apply();
		test_sym_trans_ApplyTrans();
		test_sym_trans_GetExplicitOp();
	}
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	PetscFinalize();
	return 0;
	
}