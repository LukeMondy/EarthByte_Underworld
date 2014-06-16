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

#include "private/compat/petsccompat.h"

void test_stokes_operator_K_expG( PetscTruth view  )
{
	Mat A, K,G,Gt;
	PetscInt M,N;
	Vec x, u,p;
	Vec b, f,h;
	
	
	if( view ) PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	M = 10;
	N = 8;
	
	K = MATAIJNew( PETSC_COMM_WORLD, M,M, "K" );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G" );
	
	MatFillStride( K, 1.0, 2.0 );
	MatFillStride( G, 2.0, 4.4 );
	MatTranspose( G, MAT_INITIAL_MATRIX, &Gt );
	
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes_Block( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, Gt, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	MatDestroy( & K );
	MatDestroy( & G );
	MatDestroy( & Gt );
	
	
	MatGetVecs( A, &x, &b );
	
	VecBlockGetSubVector( b, 0, &f );
	VecBlockGetSubVector( b, 1, &h );
	VecBlockRestoreSubVectors( b );
	
	VecBlockGetSubVector( x, 0, &u );
	VecBlockGetSubVector( x, 1, &p );
	VecBlockRestoreSubVectors( x );
	
	
	
	VecSet( f, 1.0 );
	VecSet( h, 2.0 );
	
	MatMult( A, b, x );
	if( view ) {
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL );
		VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	
	Stg_MatDestroy( & A );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	
}

void test_stokes_operator_K_symG( PetscTruth view )
{
	Mat A, K,G,Gt;
	PetscInt M,N,m,n;
	Vec x, u,p;
	Vec b, f,h;
	
	
	if( view ) PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
	M = 10;
	N = 8;
	
	K = MATAIJNew( PETSC_COMM_WORLD, M,M, "K" );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G" );
	
	MatGetSize( G, &M, &N );
	MatGetLocalSize( G, &m, &n );
	MatCreate( PETSC_COMM_WORLD, &Gt );
	MatSetSizes( Gt, n,m, N,M );
	MatSetType( Gt, "symtrans" );
	MatSymTransSetOperator( Gt, G );
	
	//MatCreateSymTrans( PETSC_COMM_WORLD, G, &Gt );
	
	
	MatFillStride( K, 1.0, 2.0 );
	MatFillStride( G, 2.0, 4.4 );
	
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes_Block( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	
	MatBlockSetValue( A, 0,0, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, Gt, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	Stg_MatDestroy( & K );
	Stg_MatDestroy( & G );
	Stg_MatDestroy( & Gt );
	
	
	MatGetVecs( A, &x, &b );
	
	VecBlockGetSubVector( b, 0, &f );
	VecBlockGetSubVector( b, 1, &h );
	VecBlockRestoreSubVectors( b );
	
	VecBlockGetSubVector( x, 0, &u );
	VecBlockGetSubVector( x, 1, &p );
	VecBlockRestoreSubVectors( x );
	
	
	
	VecSet( f, 1.0 );
	VecSet( h, 2.0 );
	
	MatMult( A, b, x );
	if( view ) {
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL );
		VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	
	Stg_MatDestroy( & A );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	
}





int main( int argc, char **args )
{
	int i;
	int BIG = 10000;
	PetscTruth test_mem, flg;
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	test_mem = PETSC_FALSE;
	PetscOptionsGetTruth( "", "-test_mem", &test_mem, &flg );
	
	if( test_mem ) {
		for( i=0; i<BIG; i++ ) {
			test_stokes_operator_K_symG(PETSC_FALSE);
		}
	}
	else {
		test_stokes_operator_K_expG(PETSC_TRUE);
		test_stokes_operator_K_symG(PETSC_TRUE);
	}
	
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	PetscFinalize();
	return 0;
	
}
