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
#include <petscksp.h>
#include <petscpc.h>

#include "petscext.h"
#include "petscext_vec.h"
#include "petscext_mat.h"


Mat GenerateMATSEQAIJ( MPI_Comm comm, PetscInt M, PetscInt N, const char prefix[] )
{
	Mat A;
	
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetOptionsPrefix( A, prefix );
	MatSetSizes( A, M,N, M,N );
	MatSetType( A, MATSEQAIJ );
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	return A;
}

void MatFillConstant( Mat A, PetscScalar val )
{
	PetscInt M,N, i,j;
	
	MatGetSize( A, &M, &N );
	for( i=0; i<M; i++ ) {
		for( j=0; j<N; j++ ) {
			MatSetValue( A, i,j, val, INSERT_VALUES );
		}
	}
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

void MatFillSequence( Mat A, PetscScalar start_val )
{
	PetscInt M,N, i,j;
	PetscInt cnt;
	
	MatGetSize( A, &M, &N );
	cnt = 0;
	for( i=0; i<M; i++ ) {
		for( j=0; j<N; j++ ) {
			MatSetValue( A, i,j, (PetscScalar)(start_val+cnt), INSERT_VALUES );
			cnt++;
		}
	}
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

PetscErrorCode MatBlockMergeRootBlocks( Mat B, Mat *_merge )
{
	PetscErrorCode ierr;
	PetscInt M,N, sumM, sumN;
	PetscInt MB,NB, i,j, II,JJ;
	Mat sb, A;
	MPI_Comm comm;
	PetscInt *nr, *nc;
	PetscInt *nr_offsets, *nc_offsets;
	PetscScalar val;
	PetscTruth same;
	
	sumM = 0;
	sumN = 0;
	
	
	MatGetSize( B, &MB, &NB );
	PetscMalloc( sizeof(PetscInt) * MB, &nr );
	PetscMalloc( sizeof(PetscInt) * NB, &nc );
	
	PetscMalloc( sizeof(PetscInt) * MB, &nr_offsets );
	PetscMalloc( sizeof(PetscInt) * NB, &nc_offsets );
	
	/* get max row size */
	for( i=0; i<MB; i++ ) {
		for( j=0; j<NB; j++ ) {
			sb = PETSC_NULL;
			MatBlockGetSubMatrix( B, i,j, &sb );
			if( sb != PETSC_NULL ) {
				
				PetscTypeCompare( (PetscObject)sb, "block", &same );
				if( same == PETSC_TRUE ) {
					Stg_SETERRQ( PETSC_ERR_SUP, "Operation is not recursive. Cannot merge nested blocks" );
				}
				
				MatGetSize( sb, &M, &N );
				nr[ i ]  = M;
			}
		}
	}
	
	/* get max col size */
	for( j=0; j<NB; j++ ) {
		for( i=0; i<MB; i++ ) {
			sb = PETSC_NULL;
			MatBlockGetSubMatrix( B, i,j, &sb );
			if( sb != PETSC_NULL ) {
				MatGetSize( sb, &M, &N );
				nc[ j ]  = N;
			}
		}
	}
	
	/* compute sums and offsets */
	for( i=0; i<MB; i++ ) {
		nr_offsets[i] = sumM;
		sumM = sumM + nr[i];
	}
	for( j=0; j<NB; j++ ) {
		nc_offsets[j] = sumN;
		sumN = sumN + nc[j];
	}
	
	
	PetscObjectGetComm( (PetscObject)B, &comm );
	MatCreate( comm, &A );
	MatSetSizes( A, sumM,sumN, sumM,sumN );
	MatSetType( A, MATSEQAIJ );
	
	
	/* Now build merged matrix */
	for( i=0; i<MB; i++ ) {
		for( j=0; j<NB; j++ ) {
			MatBlockGetSubMatrix( B, i,j, &sb );
			if( sb == PETSC_NULL ) continue;
			
			for( II=0; II<nr[i]; II++ ) {
				for( JJ=0; JJ<nc[j]; JJ++ ) {
					MatGetValues( sb, 1,&II, 1,&JJ, &val );
					
					MatSetValue( A, II+nr_offsets[i],JJ+nc_offsets[j], val, INSERT_VALUES );
				}
			}
		}
	}
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	
	ierr=PetscFree( nr );
	ierr=PetscFree( nc );
	
	ierr=PetscFree( nr_offsets );
	ierr=PetscFree( nc_offsets );
	
	*_merge = A;

	PetscFunctionReturn(0);
}




/*
      30    40   50 

30    [A]    [0]  [A13]
 
40    [0]    [D]  [0]

50    [A31]  [0]  [0]

*/

/*
        30+40    50 

30+40    [AA]    [BB]
 
50       [CC]    [0]


*/








PetscErrorCode test_3_field( void )
{
	Mat AA, BB, CC;
	Mat A,D , A13, A31;
	Mat S, SS; /* SS = CC * BB */
	Mat _SS;
	
	PetscPrintf( PETSC_COMM_WORLD, "\nRunning: %s \n", __func__ );
	
	/* Create sub blocks */
	A = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 3, 4, "A" );
	D = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 4, 4, "D" );
	A13 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 3, 5, "A13" );
	A31 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 5, 3, "A31" );
	
	MatFillSequence( A, 1.0 );
	MatFillSequence( D, 2.0 );
	MatFillSequence( A13, 3.0 );
	MatFillSequence( A31, 4.0 );
	
	
	/* BB */
	MatCreate( PETSC_COMM_WORLD, &BB );
	MatSetOptionsPrefix( BB, "BB" );
	MatSetSizes( BB, 2,1, 2,1 );
	MatSetType( BB, "block" );
	
	MatBlockSetValue( BB, 1,0, A13, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY);
	
	/* CC */
	MatCreate( PETSC_COMM_WORLD, &CC );
	MatSetOptionsPrefix( CC, "CC" );
	MatSetSizes( CC, 1,2, 1,2 );
	MatSetType( CC, "block" );
	
	MatBlockSetValue( CC, 0,1, A31, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(CC,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(CC,MAT_FINAL_ASSEMBLY);
	
	/* AA */
	MatCreate( PETSC_COMM_WORLD, &AA );
	MatSetOptionsPrefix( AA, "AA" );
	MatSetSizes( AA, 2,2, 2,2 );
	MatSetType( AA, "block" );
	
	MatBlockSetValue( AA, 0,0, A, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( AA, 1,1, D, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY);
	
	{
		Mat AA_seqaij;
		MatBlockMergeRootBlocks( AA, &AA_seqaij );
		PetscPrintf( PETSC_COMM_WORLD, "Merged version of AA \n");
		PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
		MatView( AA_seqaij, PETSC_VIEWER_STDOUT_WORLD );
	}
	
	
	
	/* Do std multiplication */
	PetscPrintf( PETSC_COMM_WORLD, "Std mat mat mult \n");
	MatMatMult( A31, A13, MAT_INITIAL_MATRIX, 2.0, &S );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( S, PETSC_VIEWER_STDOUT_WORLD );
	
	
	/* Do block multiplication */
	PetscPrintf( PETSC_COMM_WORLD, "Block mat mat mult \n");
	MatMatMult( CC, BB, MAT_INITIAL_MATRIX, 2.0, &SS );
//	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
//	MatView( SS, PETSC_VIEWER_STDOUT_WORLD );
	
	MatBlockGetSubMatrix( SS, 0,0, &_SS );
	MatBlockRestoreSubMatrices( SS );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( _SS, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_MatDestroy( SS );
	Stg_MatDestroy( A );
	Stg_MatDestroy( D );
	Stg_MatDestroy( A13 );
	Stg_MatDestroy( A31 );
	
	PetscFunctionReturn(0);	
}







PetscErrorCode test_2_block( void )
{
	PetscErrorCode ierr;
	Mat A13,A23 , A31,A32;
	Mat G, D;
	Mat Gmerged, Dmerged;
	Mat DG;
	Mat DmGm, _DG;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\nRunning: %s \n", __func__ );
	
	/* Create sub blocks */
	A13 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 3, 5, "A13" );
	A23 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 4, 5, "A23" );
	MatFillSequence( A13, 3.0 );
	MatFillSequence( A23, 4.0 );
	
	A31 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 5, 3, "A31" );
	A32 = GenerateMATSEQAIJ( PETSC_COMM_WORLD, 5, 4, "A32" );
	MatFillSequence( A31, 5.0 );
	MatFillSequence( A32, 6.0 );
	
	
	
	/* G */
	MatCreate( PETSC_COMM_WORLD, &G );
	MatSetOptionsPrefix( G, "G" );
	MatSetSizes( G, 2,1, 2,1 );
	MatSetType( G, "block" );
	
	ierr=MatBlockSetValue( G, 0,0, A13, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);CHKERRQ(ierr);
	ierr=MatBlockSetValue( G, 1,0, A23, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);CHKERRQ(ierr);
	
	MatAssemblyBegin(G,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(G,MAT_FINAL_ASSEMBLY);
	
	/* D */
	MatCreate( PETSC_COMM_WORLD, &D );
	MatSetOptionsPrefix( D, "D" );
	MatSetSizes( D, 1,2, 1,2 );
	MatSetType( D, "block" );
	
	MatBlockSetValue( D, 0,0, A31, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	MatBlockSetValue( D, 0,1, A32, DIFFERENT_NONZERO_PATTERN, INSERT_VALUES);
	
	MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY);
	
	MatBlockMergeRootBlocks( G, &Gmerged );
	MatBlockMergeRootBlocks( D, &Dmerged );
	
	/* Do std multiplication */
	PetscPrintf( PETSC_COMM_WORLD, "Std mat mat mult \n");
	MatMatMult( Dmerged, Gmerged, MAT_INITIAL_MATRIX, 2.0, &DmGm );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( DmGm, PETSC_VIEWER_STDOUT_WORLD );
	
	/* Do block multiplication */
	PetscPrintf( PETSC_COMM_WORLD, "Block mat mat mult \n");
	MatMatMult( D, G, MAT_INITIAL_MATRIX, 2.0, &DG );
	MatSetOptionsPrefix( DG, "DG" );
	
	MatBlockGetSubMatrix( DG, 0,0, &_DG );
	MatBlockRestoreSubMatrices( DG );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( _DG, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_MatDestroy( D );
	Stg_MatDestroy( G );
	Stg_MatDestroy( Dmerged );
	Stg_MatDestroy( Gmerged );
	
	Stg_MatDestroy( DmGm );
	Stg_MatDestroy( DG );

	PetscFunctionReturn(0);	
}




int main( int argc, char **args )
{
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtVecRegisterAll();
	PetscExtMatRegisterAll();
	
	
	
	test_2_block();
	test_3_field();
	
	PetscExtVecRegisterDestroyAll();
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
