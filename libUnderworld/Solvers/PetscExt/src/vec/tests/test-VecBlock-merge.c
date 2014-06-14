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

#include "petscext_vec.h"





/*
X = ( a b )^T


*/

void test_merge_seq( void )
{
	Vec X, a,b, mX;
	PetscInt index;
	PetscReal val;
	
	
	PetscPrintf( PETSC_COMM_WORLD, "\n\n============== %s ==============\n", __func__ );
	
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2, 2 );
	VecSetType( X, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &a );
	VecSetSizes( a, PETSC_DECIDE, 2 );
	VecSetType( a, VECSEQ );
	VecSet( a, 1.0 );
	
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, PETSC_DECIDE, 3 );
	VecSetType( b, VECSEQ );
	VecSet( b, 2.0 );
	
	
	/* assemble X */
	VecBlockSetValue( X, 0, a, INSERT_VALUES );
	VecBlockSetValue( X, 1, b, INSERT_VALUES );
	Stg_VecDestroy( & a );
	Stg_VecDestroy( & b );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	
//	PetscPrintf( PETSC_COMM_WORLD, "X \n");
//	VecView( X, PETSC_VIEWER_STDOUT_WORLD );
	
	
	mX = PETSC_NULL;
	VecBlockMergeSubVectors( X, INSERT_VALUES, VECSEQ, &mX );
	VecBlockMergeSubVectors( X, ADD_VALUES, VECSEQ, &mX );

	PetscPrintf( PETSC_COMM_WORLD, "mX \n");
	VecView( mX, PETSC_VIEWER_STDOUT_WORLD );
	
	Stg_VecDestroy( & X );
	Stg_VecDestroy( & mX );
	
}

void test_merge( int rank )
{
	Vec X, a,b, mX;
	PetscInt index;
	PetscReal val;
	
	
	
	if( rank ==0 ) PetscPrintf( PETSC_COMM_WORLD, "\n\n============== %s ==============\n", __func__ );
	
	
	VecCreate( PETSC_COMM_WORLD, &X );
	VecSetSizes( X, 2, 2 );
	VecSetType( X, "block" );
	
	VecCreate( PETSC_COMM_WORLD, &a );
	VecSetSizes( a, PETSC_DECIDE, 5 );
	VecSetType( a, VECMPI );
	VecSet( a, 1.0 );
	
	VecCreate( PETSC_COMM_WORLD, &b );
	VecSetSizes( b, PETSC_DECIDE, 6 );
	VecSetType( b, VECMPI );
	VecSet( b, 2.0 );
	
	
	/* assemble X */
	VecBlockSetValue( X, 0, a, INSERT_VALUES );
	VecBlockSetValue( X, 1, b, INSERT_VALUES );
	Stg_VecDestroy( & a );
	Stg_VecDestroy( & b );
	VecAssemblyBegin(X);
	VecAssemblyEnd(X);
	
	//if( rank ==0 ) PetscPrintf( PETSC_COMM_WORLD, "X \n");
	//VecView( X, PETSC_VIEWER_STDOUT_WORLD );
	

	mX = PETSC_NULL;
	VecBlockMergeSubVectors( X, INSERT_VALUES, VECMPI, &mX );

	if( rank ==0 ) PetscPrintf( PETSC_COMM_WORLD, "mX \n");
	VecView( mX, PETSC_VIEWER_STDOUT_WORLD );
	
	Stg_VecDestroy( & mX );
	Stg_VecDestroy( & X );
}




int main( int argc, char **args )
{
	int i;
	int BIG;
	int nproc, rank;
	
	PetscInitialize( &argc, &args,(char *)0, PETSC_NULL );
	PetscExtVecRegisterAll();
	
	MPI_Comm_size( PETSC_COMM_WORLD, &nproc );
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	
	
	for( i=0; i<1; i++ ) {
		
		if( nproc == 1 ) {
			test_merge_seq();
		}
		else {
			test_merge(rank);
		}
	}
	
	PetscExtVecRegisterDestroyAll();
	PetscFinalize();
	return 0;
}
