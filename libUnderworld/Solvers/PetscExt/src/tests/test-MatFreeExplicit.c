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
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "petscext_helpers.h"
#include "petscext_utils.h"
#include "petscext_mat.h"



void test_MatFreeExplicit( void )
{
	Mat A,AE;
	PetscInt M,N;

	M = 9;
	N = 4;
	A = MATAIJNew( PETSC_COMM_WORLD, M,N, "A_" );
	MatFillStride( A, 1.0, 2.2 );

	MatViewContents(A,"A_");

	MatFreeComputeExplicitOperator( A, MATAIJ, 1.0e-6, &AE );
	MatViewContents(AE,"AE_");

	MatDestroy( &&A);
	MatDestroy( &&AE);
}

void test_MatFreeSchurExplicit_11( void )
{
	Mat K,G,D,C,schur,S;
	PetscInt M,N;
	KSP ksp;
	PC pc;

	M = 9;
	N = 5;

	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G_" );
	MatFillStride( G, 1.0, 2.2 );

	D = MATAIJNew( PETSC_COMM_WORLD, N,M, "D_" );
	MatFillStride( D, 2.0, 1.2 );

	C = MATAIJNew( PETSC_COMM_WORLD, N,N, "C_" );
	MatFillStride( C, 1.0, 0.0 );
	
	K = MATAIJIdentityNew( PETSC_COMM_WORLD, M );

        MatCreateSchur( PETSC_COMM_WORLD, K,G,D,C, PETSC_NULL, "MatSchur_A11", &schur );
        MatDestroy( & K );        MatDestroy( & G );
	MatDestroy( & D );        MatDestroy( & C );

        MatSchurGetKSP( schur, &ksp );
        KSPSetType( ksp, "preonly" );
        KSPGetPC( ksp, &pc );
        PCSetType( pc, "jacobi" );

        MatAssemblyBegin( schur, MAT_FINAL_ASSEMBLY );
        MatAssemblyEnd( schur, MAT_FINAL_ASSEMBLY );
	
	MatFreeComputeExplicitOperator( schur, MATAIJ, 1.0e-6, &S );
	MatViewContents(S,"S_");

	MatDestroy( &schur);
	MatDestroy( &S);
}

void test_MatFreeSchurExplicit_22( void )
{
	Mat K,G,D,C,schur,S;
	PetscInt M,N;
	KSP ksp;
	PC pc;

	M = 9;
	N = 5;

	G = MATAIJNew( PETSC_COMM_WORLD, M,N, "G_" );
	MatFillStride( G, 1.0, 2.2 );

	D = MATAIJNew( PETSC_COMM_WORLD, N,M, "D_" );
	MatFillStride( D, 2.0, 1.2 );

	K = MATAIJNew( PETSC_COMM_WORLD, M,M, "K_" );
	MatFillStride( K, 3.0, 5.76 );
	
	C = MATAIJIdentityNew( PETSC_COMM_WORLD, N );


        MatCreateSchur( PETSC_COMM_WORLD, K,G,D,C, PETSC_NULL, "MatSchur_A22", &schur );
        MatDestroy( & K );        MatDestroy( & G );
	MatDestroy( & D );        MatDestroy( & C );

        MatSchurGetKSP( schur, &ksp );
        KSPSetType( ksp, "preonly" );
        KSPGetPC( ksp, &pc );
        PCSetType( pc, "jacobi" );

        MatAssemblyBegin( schur, MAT_FINAL_ASSEMBLY );
        MatAssemblyEnd( schur, MAT_FINAL_ASSEMBLY );
	
	MatFreeComputeExplicitOperator( schur, MATAIJ, 1.0e-6, &S );
	MatViewContents(S,"S_");

	MatDestroy( &schur);
	MatDestroy( &S);
}

int main( int argc, char **args )
{ 
        PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	PetscExtMatRegisterAll();
        
        test_MatFreeExplicit();
        test_MatFreeSchurExplicit_11();
        test_MatFreeSchurExplicit_22();
        
	PetscExtMatRegisterDestroyAll();
        PetscFinalize();
        return 0;        
}

