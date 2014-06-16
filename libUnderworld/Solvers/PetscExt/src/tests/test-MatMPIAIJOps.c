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

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "petscext_utils.h"


void test_MatSeqAIJOps( void )
{
	Mat A;
	Vec y;
	PetscInt *idx,m,n,i;
	PetscScalar *_y;
	PetscMPIInt rank;
	
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	if (rank!=0) { return; }
	

	MatAIJLoad_MatrixMarket( PETSC_COMM_SELF, "input/e05r0000.mtx", &A );
	
	MatGetVecs( A, PETSC_NULL, &y );
	
	/*  */
	MatGetLocalSize( A, &m,&n );
	PetscMalloc( sizeof(PetscInt)*m, &idx );
	
	MatGetRowMax(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMax \n" );
	if (rank==0) {		VecView( y, PETSC_VIEWER_STDOUT_SELF );		}
	/*
	VecGetArray( y, &_y );
	for (i=0; i<m; i++) {
		PetscPrintf( PETSC_COMM_WORLD, "  [%d]: %f \n", i, PetscRealPart(_y[i]) );
	}
	VecRestoreArray( y, &_y );
	*/
	
	MatGetRowMin(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMin \n" );
	if (rank==0) {		VecView( y, PETSC_VIEWER_STDOUT_SELF );		}
	/*
	VecGetArray( y, &_y );
	for (i=0; i<m; i++) {
		PetscPrintf( PETSC_COMM_WORLD, "  [%d]: %f \n", i, PetscRealPart(_y[i]) );
	}
	VecRestoreArray( y, &_y );
	*/
	
	MatGetRowMaxAbs(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMaxAbs \n" );
	if (rank==0) {		VecView( y, PETSC_VIEWER_STDOUT_SELF );		}
	/*
	VecGetArray( y, &_y );
	for (i=0; i<m; i++) {
		PetscPrintf( PETSC_COMM_WORLD, "  [%d]: %f \n", i, PetscRealPart(_y[i]) );
	}
	VecRestoreArray( y, &_y );
	*/
	
	PetscFree( idx );
	MatDestroy( & A );
	Stg_VecDestroy( & y );
}


void test_MatMPIAIJOps( void )
{
	Mat A;
	Vec y;
	PetscInt *idx,m,n;
	
	
	MatAIJLoad_MatrixMarket( PETSC_COMM_WORLD, "input/e05r0000.mtx", &A );
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	MatGetVecs( A, PETSC_NULL, &y );
	
	/*  */
	MatGetLocalSize( A, &m,&n );
	PetscMalloc( sizeof(PetscInt)*m, &idx );
	
	PetscExtMatGetRowMax_MPIAIJ(A,y,PETSC_NULL);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMax \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscExtMatGetRowMin_MPIAIJ(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMin \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	PetscExtMatGetRowMaxAbs_MPIAIJ(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMaxAbs \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	
	PetscFree( idx );
	MatDestroy( & A );
	Stg_VecDestroy( & y );
}


void test_MatMPIAIJOps2( void )
{
	Mat A;
	Vec y;
	PetscInt *idx,m,n;
	
	
	MatAIJLoad_MatrixMarket( PETSC_COMM_WORLD, "input/e05r0000.mtx", &A );
	PetscExtLoadMPIAIJOperations(A);
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	MatGetVecs( A, PETSC_NULL, &y );
	
	/*  */
	MatGetLocalSize( A, &m,&n );
	PetscMalloc( sizeof(PetscInt)*m, &idx );
	
	MatGetRowMax(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMax \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	MatGetRowMin(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMin \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	MatGetRowMaxAbs(A,y,idx);
	PetscPrintf( PETSC_COMM_WORLD, "MatGetRowMaxAbs \n" );
	VecView( y, PETSC_VIEWER_STDOUT_WORLD );
	
	
	PetscFree( idx );
	MatDestroy( & A );
	Stg_VecDestroy( & y );
}





int main( int argc, char **args )
{
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	
	test_MatSeqAIJOps();
//	test_MatMPIAIJOps();
	test_MatMPIAIJOps2();
	
	
	
	PetscFinalize();
	return 0;
	
}
