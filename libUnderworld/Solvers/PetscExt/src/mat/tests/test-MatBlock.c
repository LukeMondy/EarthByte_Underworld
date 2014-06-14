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


void test_solve( void )
{
	Mat A11, A12,A21, A;
	KSP ksp;
	PC pc;
	Vec b,x , f,h, diag, x1,x2;
	int n, np, i,j;
	const MatType type;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "%s \n", __func__ );
	
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
	
	// A12
	MatCreate( PETSC_COMM_WORLD, &A12 );
	MatSetSizes( A12, PETSC_DECIDE, PETSC_DECIDE, n, np );
	MatSetType( A12, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A12, np, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A12, i,j, (double)(i+j*n), INSERT_VALUES );
		}
	}
	MatSetValue( A12, 2,1, (double)(4), INSERT_VALUES );
	MatAssemblyBegin( A12, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A12, MAT_FINAL_ASSEMBLY);
	
	// A21
	MatTranspose( A12, MAT_INITIAL_MATRIX, &A21 );
	
	
	
	/* Create block matrix */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, 2,2, 2,2 );
	MatSetType( A, "block" );
	MatGetType( A, &type );
	
	
	MatBlockSetValue( A, 0,0, A11, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 0,1, A12, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES );
	MatBlockSetValue( A, 1,0, A21, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	

        {
                PetscInt M,N,m,n;
                PetscInt i,s,e;
                const PetscInt *r;
                PetscMPIInt np,rank;
        

                MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
                MPI_Comm_size( PETSC_COMM_WORLD, &np );

                MatGetSize( A, &M,&N);
        //        printf("M,N = %d,%d \n", M,N);

                MatGetLocalSize( A, &m,&n );
        //        printf("m,n = %d,%d \n", m,n );

                MatGetOwnershipRange( A, &s,&e);
        //        printf("s,e = %d,%d \n", s,e );

        //      MatGetOwnershipRanges( S, &r );
        //      for( i=0; i<np+1; i++ ) {
        //              printf("r[%d]=%d \n", i, r[i] );
        //      }
        }


	
	/* Create vectors */
	MatGetVecs( A12, &h, &f );
	
	VecSet( f, 1.0 );
	VecSet( h, 0.0 );
	
	
	/* Create block vector */
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, 2, 2 );
	VecSetType( b, "block" );
	
	VecBlockSetValue( b, 0, f, INSERT_VALUES );
	VecBlockSetValue( b, 1, h, INSERT_VALUES );
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	
	VecDuplicate( b, &x );
	
	
	KSPCreate( PETSC_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPSetType( ksp, "gmres" );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, "none" );
	KSPSetFromOptions( ksp );
	
	KSPSolve( ksp, b, x );
	
	
	VecBlockGetSubVector( x, 0, &x1 );
	VecBlockGetSubVector( x, 1, &x2 );
	VecBlockRestoreSubVectors( x );
	
	PetscPrintf( PETSC_COMM_WORLD, "x1 \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
	VecView( x1, PETSC_VIEWER_STDOUT_WORLD );
	PetscPrintf( PETSC_COMM_WORLD, "x2 \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
	VecView( x2, PETSC_VIEWER_STDOUT_WORLD );
	
	
	/*
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO_DETAIL );
	VecView( x, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	KSPDestroy( & ksp );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & b );
	Stg_MatDestroy( A11 );
	Stg_MatDestroy( A12 );
	Stg_MatDestroy( A21 );
	Stg_VecDestroy( & f );
	Stg_VecDestroy( & h );
	
	Stg_MatDestroy( A );
	
}

void test_3_field( void )
{
	Mat THREE_FIELD_OPERATOR;
	Mat AA, BB, CC;
	Mat A,B,C , Bp, Cp;
	
	
	/* BB */
	MatCreate( PETSC_COMM_WORLD, &Bp );
//	PetscObjectSetName( (PetscObject)Bp, "Bp" );
	MatSetOptionsPrefix( Bp, "Bp" );
	MatSetSizes( Bp, 40,50, 40,50 );
	MatSetType( Bp, MATSEQAIJ );
	MatAssemblyBegin(Bp,MAT_FINAL_ASSEMBLY);	MatAssemblyEnd(Bp,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "Bp \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( Bp, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &BB );
	MatSetOptionsPrefix( BB, "BB" );
	MatSetSizes( BB, 2,1, 2,1 );
	MatSetType( BB, "block" );
	MatBlockSetValue( BB, 1,0, Bp, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "BB \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( BB, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	/* CC */
	MatCreate( PETSC_COMM_WORLD, &Cp );
	MatSetOptionsPrefix( Cp, "Cp" );
	MatSetSizes( Cp, 50,40, 50,40 );
	MatSetType( Cp, MATSEQAIJ );
	MatAssemblyBegin(Cp,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Cp,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "Cp \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( Cp, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &CC );
	MatSetOptionsPrefix( CC, "CC" );
	MatSetSizes( CC, 1,2, 1,2 );
	MatSetType( CC, "block" );
	MatBlockSetValue( CC, 0,1, Cp, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(CC,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(CC,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "CC \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( CC, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	/* AA */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetOptionsPrefix( A, "A" );
	MatSetSizes( A, 30,30, 30,30 );
	MatSetType( A, MATSEQAIJ );
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "A \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &B );
	MatSetOptionsPrefix( B, "B" );
	MatSetSizes( B, 30,40, 30,40 );
	MatSetType( B, MATSEQAIJ );
	MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "B \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( B, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &C );
	MatSetOptionsPrefix( C, "C" );
	MatSetSizes( C, 40,30, 40,30 );
	MatSetType( C, MATSEQAIJ );
	MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "C \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( C, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	MatCreate( PETSC_COMM_WORLD, &AA );
	MatSetOptionsPrefix( AA, "AA" );
	MatSetSizes( AA, 2,2, 2,2 );
	MatSetType( AA, "block" );
	MatBlockSetValue( AA, 0,0, A, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( AA, 0,1, B, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( AA, 1,0, C, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "AA \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( AA, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	MatCreate( PETSC_COMM_WORLD, &THREE_FIELD_OPERATOR );
	MatSetOptionsPrefix( THREE_FIELD_OPERATOR, "3-field" );
//	PetscObjectSetName( (PetscObject)THREE_FIELD_OPERATOR, "3-field" );
	MatSetSizes( THREE_FIELD_OPERATOR, 2,2, 2,2 );
	MatSetType( THREE_FIELD_OPERATOR, "block" );
	MatBlockSetValue( THREE_FIELD_OPERATOR, 0,0, AA, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( THREE_FIELD_OPERATOR, 0,1, BB, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( THREE_FIELD_OPERATOR, 1,0, CC, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(THREE_FIELD_OPERATOR,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(THREE_FIELD_OPERATOR,MAT_FINAL_ASSEMBLY);
	
	
	
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( THREE_FIELD_OPERATOR, PETSC_VIEWER_STDOUT_WORLD );
	
	
}


/*
      30    40   50 

30    [A]   [B]  [0]

40    [C]   [0]  [Bp]

50    [D]   [Cp] [0]

*/


void test_mat_get_vecs( void )
{
	Mat THREE_FIELD_OPERATOR;
	Mat AA, BB, CC;
	Mat A,B,C , Bp, Cp,D;
	Vec AAL, AAR, TFL, TFR;
	Vec cl,cr;
	PetscInt s;
	
	/* BB */
	MatCreate( PETSC_COMM_WORLD, &Bp );
//	PetscObjectSetName( (PetscObject)Bp, "Bp" );
	MatSetOptionsPrefix( Bp, "Bp" );
	MatSetSizes( Bp, 40,50, 40,50 );
	MatSetType( Bp, MATSEQAIJ );
	MatAssemblyBegin(Bp,MAT_FINAL_ASSEMBLY);	MatAssemblyEnd(Bp,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "Bp \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( Bp, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &BB );
	MatSetOptionsPrefix( BB, "BB" );
	MatSetSizes( BB, 2,1, 2,1 );
	MatSetType( BB, "block" );
	MatBlockSetValue( BB, 1,0, Bp, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "BB \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( BB, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	/* CC */
	MatCreate( PETSC_COMM_WORLD, &Cp );
	MatSetOptionsPrefix( Cp, "Cp" );
	MatSetSizes( Cp, 50,40, 50,40 );
	MatSetType( Cp, MATSEQAIJ );
	MatAssemblyBegin(Cp,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Cp,MAT_FINAL_ASSEMBLY);
	
	MatCreate( PETSC_COMM_WORLD, &D );
	MatSetOptionsPrefix( D, "D" );
	MatSetSizes( D, 50,30, 50,30 );
	MatSetType( D, MATSEQAIJ );
	MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY);
	
	
	/*
	MatGetVecs( Cp, &cr, &cl );
	VecGetSize( cr, &s );
	PetscPrintf( PETSC_COMM_WORLD, "cr [%d] \n", s);
	VecGetSize( cl, &s );
	PetscPrintf( PETSC_COMM_WORLD, "cl [%d] \n", s);
	*/
	
	
	
	/*	
	PetscPrintf( PETSC_COMM_WORLD, "Cp \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( Cp, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &CC );
	MatSetOptionsPrefix( CC, "CC" );
	MatSetSizes( CC, 1,2, 1,2 );
	MatSetType( CC, "block" );
	MatBlockSetValue( CC, 0,0, D, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( CC, 0,1, Cp, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(CC,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(CC,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "CC \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( CC, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	/* cl \in 50 x 1, cr \in { 30 x 1, 40 x 1 } */
	MatGetVecs( CC, &cr, &cl ); /* cl = CC cr */
	PetscPrintf( PETSC_COMM_WORLD, "cr \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( cr, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscPrintf( PETSC_COMM_WORLD, "cl \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( cl, PETSC_VIEWER_STDOUT_WORLD );
	
	
	
	
	/* AA */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetOptionsPrefix( A, "A" );
	MatSetSizes( A, 30,30, 30,30 );
	MatSetType( A, MATSEQAIJ );
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "A \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &B );
	MatSetOptionsPrefix( B, "B" );
	MatSetSizes( B, 30,40, 30,40 );
	MatSetType( B, MATSEQAIJ );
	MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "B \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( B, PETSC_VIEWER_STDOUT_WORLD );
*/	
	MatCreate( PETSC_COMM_WORLD, &C );
	MatSetOptionsPrefix( C, "C" );
	MatSetSizes( C, 40,30, 40,30 );
	MatSetType( C, MATSEQAIJ );
	MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "C \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( C, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	MatCreate( PETSC_COMM_WORLD, &AA );
	MatSetOptionsPrefix( AA, "AA" );
	MatSetSizes( AA, 2,2, 2,2 );
	MatSetType( AA, "block" );
	MatBlockSetValue( AA, 0,0, A, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( AA, 0,1, B, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( AA, 1,0, C, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY);
/*	
	PetscPrintf( PETSC_COMM_WORLD, "AA \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( AA, PETSC_VIEWER_STDOUT_WORLD );
*/	
	
	MatCreate( PETSC_COMM_WORLD, &THREE_FIELD_OPERATOR );
	MatSetOptionsPrefix( THREE_FIELD_OPERATOR, "3-field" );
//	PetscObjectSetName( (PetscObject)THREE_FIELD_OPERATOR, "3-field" );
	MatSetSizes( THREE_FIELD_OPERATOR, 2,2, 2,2 );
	MatSetType( THREE_FIELD_OPERATOR, "block" );
	MatBlockSetValue( THREE_FIELD_OPERATOR, 0,0, AA, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( THREE_FIELD_OPERATOR, 0,1, BB, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( THREE_FIELD_OPERATOR, 1,0, CC, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatAssemblyBegin(THREE_FIELD_OPERATOR,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(THREE_FIELD_OPERATOR,MAT_FINAL_ASSEMBLY);
	
	
	
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( THREE_FIELD_OPERATOR, PETSC_VIEWER_STDOUT_WORLD );
	
	
	PetscPrintf( PETSC_COMM_WORLD, "***********************************\n");
	MatGetVecs( AA, &AAR,&AAL );
	
	PetscPrintf( PETSC_COMM_WORLD, "AAR \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( AAR, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscPrintf( PETSC_COMM_WORLD, "AAL \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( AAL, PETSC_VIEWER_STDOUT_WORLD );
	
	MatGetVecs( THREE_FIELD_OPERATOR, &TFR,&TFL );

	PetscPrintf( PETSC_COMM_WORLD, "TFR \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( TFR, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscPrintf( PETSC_COMM_WORLD, "TFL \n");
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	VecView( TFL, PETSC_VIEWER_STDOUT_WORLD );
	
}

int main( int argc, char **args )
{
	int i;
	int BIG=1;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	
	test_solve();
//	test_3_field();
//	test_mat_get_vecs();
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
