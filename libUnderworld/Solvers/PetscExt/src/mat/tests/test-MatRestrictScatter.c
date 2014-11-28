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


void test_MatRestrictScatter( void )
{
	Mat A,B;
	IS row,col;
	PetscInt n,np, i,j, s,e;
	Vec L,R;
	
	n = 13;
	np = 10;
	
	// A
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, n, np );
	MatSetType( A, MATAIJ );
	MatSeqAIJSetPreallocation( A, np, PETSC_NULL );
	
	for( i=0; i<n; i++ ) {
		for( j=0; j<np; j++ ) {
			MatSetValue( A, i,j, (double)(i+j*n+1.0), INSERT_VALUES );
		}
	}
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY);
	
	/* IS */
	MatGetOwnershipRange( A, &s,&e );
	ISCreateStride( PETSC_COMM_WORLD, 2, s, 1, &row );
	ISCreateStride( PETSC_COMM_WORLD, 2, s, 1, &col );
	
	MatCreateRestrictScatter( PETSC_COMM_WORLD, A, row,col, &B );
	MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_DENSE );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	
	MatGetVecs( B, &R, &L );
	VecSet(R,1.0);
	MatMult( B,R, L );
	
	PetscPrintf( PETSC_COMM_WORLD, "L = \n");
	VecView( L, PETSC_VIEWER_STDOUT_WORLD );
	
	
	Stg_VecDestroy( &L);
	Stg_VecDestroy( &R);
	Stg_ISDestroy(row);
	Stg_ISDestroy(col);
	Stg_MatDestroy( A );
	Stg_MatDestroy( B );
}



int main( int argc, char **args )
{
	int i;
	int BIG=1;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	
	PetscExtMatRegisterAll();
	
	for( i=0; i<BIG; i++ ) {
		test_MatRestrictScatter();
	}
	
	PetscExtMatRegisterDestroyAll();
	
	PetscFinalize();
	return 0;
}
