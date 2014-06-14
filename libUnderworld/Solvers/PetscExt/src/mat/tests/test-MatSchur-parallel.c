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

#include "stdio.h"
#include "stdlib.h"

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>

#include "petscext_helpers.h"
#include "petscext_utils.h"
#include "petscext_vec.h"
#include "petscext_mat.h"


PetscErrorCode _MatView(Mat A,const char name[])
{
	PetscPrintf( PETSC_COMM_WORLD, "mat: %s\n", name );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	PetscFunctionReturn(0);
}
PetscErrorCode _VecView(Vec x,const char name[])
{
	PetscPrintf( PETSC_COMM_WORLD, "vec: %s\n", name );
	VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	PetscFunctionReturn(0);
}

/*
Build a non symmetric system with A22. Check MatMult, MatMultTranspose
*/
void test_schur_11_test2( void )
{
	Mat A_bk;
	Mat K,G,D,C;
	PetscInt M,N;
	Mat schur;
	Vec p, fhat;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "Test: %s \n", __func__ );
	

	M = 7;
	N = 5;
	
	K = MATAIJIdentityNew( PETSC_COMM_WORLD, M );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G_" );
	D = MATAIJNew( PETSC_COMM_WORLD, N,M, "D_" );
	C = MATAIJNew( PETSC_COMM_WORLD, N,N, "C_" );
	
	MatFillStride( G, 2.0, 4.4 );
	MatFillStride( D, 5.0, 1.3 );
	MatFillStride( C, 1.0, 1.0 );
	
	
	MatCreate( PETSC_COMM_WORLD, &A_bk );
	MatSetSizes( A_bk, 2,2, 2,2 );
	MatSetType( A_bk, "block" );
	
	MatBlockSetValue( A_bk, 0,0, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A_bk, 0,1, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A_bk, 1,0, D,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A_bk, 1,1, C,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A_bk,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_bk,MAT_FINAL_ASSEMBLY);
	Stg_MatDestroy( K );	Stg_MatDestroy( G );
	Stg_MatDestroy( D );	Stg_MatDestroy( C );
	
	
	MatCreateSchurFromBlock( A_bk, PETSC_NULL, "MatSchur_A11", &schur );
	MatGetVecs( schur, &p, &fhat );
	
	VecSet( p, 1.0 );
	
	MatAssemblyBegin( schur, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( schur, MAT_FINAL_ASSEMBLY );
	
	MatMult( schur, p, fhat );
	_VecView(fhat,"S{1}");
	MatMultTranspose( schur, p, fhat );
	_VecView(fhat,"S^T{1}");

	
	
	Stg_MatDestroy( schur );
	
	Stg_VecDestroy( & p );
	Stg_VecDestroy( & fhat );
	Stg_MatDestroy( A_bk );
}

/*
Build a non symmetric system with A11 and A22. 
Set scale.
Set pc type lu
Check MatMult, MatMultTranspose
*/
void test_schur_22_test3( void )
{
	Mat A_bk;
	Mat K,G,D,C;
	PetscInt M,N;
	Mat schur;
	Vec u, f;
	KSP ksp;
	PC pc;
	
	PetscPrintf( PETSC_COMM_WORLD, "Test: %s \n", __func__ );

	M = 7;
	N = 5;
	
	K = MATAIJIdentityNew( PETSC_COMM_WORLD, M );
	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G_" );
	D = MATAIJNew( PETSC_COMM_WORLD, N,M, "D_" );
	C = MATAIJIdentityNew( PETSC_COMM_WORLD, N );
	
	MatFillStride( G, 2.0, 4.4 );
	MatFillStride( D, 5.0, 1.3 );
	MatScale( C, 2.0 );
	
	
	MatCreate( PETSC_COMM_WORLD, &A_bk );
	MatSetSizes( A_bk, 2,2, 2,2 );
	MatSetType( A_bk, "block" );
	
	MatBlockSetValue( A_bk, 0,0, K,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A_bk, 0,1, G,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A_bk, 1,0, D,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( A_bk, 1,1, C,  DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A_bk,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_bk,MAT_FINAL_ASSEMBLY);
	Stg_MatDestroy( K );	Stg_MatDestroy( G );
	Stg_MatDestroy( D );	Stg_MatDestroy( C );
	
	
	MatCreateSchurFromBlock( A_bk, 4.9, "MatSchur_A22", &schur );
	MatSchurGetKSP( schur, &ksp );
	
	MatAssemblyBegin( schur, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( schur, MAT_FINAL_ASSEMBLY );
	
	MatGetVecs( schur, &u, &f );
	VecSet( u, 1.0 );
	
	MatMult( schur, u, f );
	_VecView(f,"S{1}");
	MatMultTranspose( schur, u, f );
	_VecView(f,"S^T{1}");
	

	{
		PetscInt M,N,m,n;
		PetscInt i,s,e;
		const PetscInt *r;
		PetscMPIInt np,rank;
		Mat S;

		S = schur;
	//	MatCreateNormal(K,&S);
	//	MatView( S, PETSC_VIEWER_STDOUT_WORLD );

		MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
		MPI_Comm_size( PETSC_COMM_WORLD, &np );

		MatGetSize( S, &M,&N);
	//	printf("M,N = %d,%d \n", M,N);

		MatGetLocalSize( S, &m,&n );
	//	printf("m,n = %d,%d \n", m,n );

		MatGetOwnershipRange( S, &s,&e);
	//	printf("s,e = %d,%d \n", s,e );

	//	MatGetOwnershipRanges( S, &r );
	//	for( i=0; i<np+1; i++ ) {
	//		printf("r[%d]=%d \n", i, r[i] );
	//	}
	}
	
	
	Stg_MatDestroy( schur );	
	Stg_VecDestroy( & u );
	Stg_VecDestroy( & f );
	Stg_MatDestroy( A_bk );
}



int main(int argc,char **argv)
{
	int i;
	Vec x;
	
	PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();

	
	for( i=0; i<1; i++ ) {
		test_schur_11_test2();PetscPrintf(PETSC_COMM_WORLD,"\n");
		test_schur_22_test3();PetscPrintf(PETSC_COMM_WORLD,"\n");
	}	

        
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	PetscFinalize();
	return 0;
}
