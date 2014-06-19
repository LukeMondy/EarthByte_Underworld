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

#include "petscext_helpers.h"
#include "petscext_utils.h"


void test_VecLoadMatrixMarket( void )
{
	Vec x,y;
	PetscReal v_norm;
	PetscViewer viewer;
	
	
	VecLoad_MatrixMarket( PETSC_COMM_WORLD, "input/e05r0000_rhs1.mtx", PETSC_NULL, &x );
	VecNorm( x, NORM_2, &v_norm );
	PetscPrintf( PETSC_COMM_WORLD, "norm({x}) = %5.5e \n\n", v_norm );
	
	PetscPrintf( PETSC_COMM_WORLD, "Dump as MM \n");
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, "petsc_rhs1.mtx", &viewer );
	VecView_MatrixMarket( x, viewer );
	Stg_PetscViewerDestroy( & viewer );
	
	Stg_VecDestroy( & x );
	
	VecLoad_MatrixMarket( PETSC_COMM_WORLD, "petsc_rhs1.mtx", PETSC_NULL, &y );
	VecNorm( y, NORM_2, &v_norm );
	PetscPrintf( PETSC_COMM_WORLD, "norm({y}) = %5.5e \n\n", v_norm );
	
	Stg_VecDestroy( & y );
}

void test_MatAIJLoadMatrixMarket( void )
{
	Mat A, An;
	PetscViewer viewer;
	Vec x,y;
	PetscScalar v_norm, m_norm;
	
	__MatLoadMatrixMarket( PETSC_COMM_WORLD, "input/e05r0000.mtx", MATAIJ, &A );
	
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO );
	MatView( A, PETSC_VIEWER_STDOUT_WORLD );
	
	
	MatNorm( A, NORM_1, &m_norm );
	
	MatGetVecs( A, &x, &y );
	VecSet( x, 1.0 );
	MatMult( A,x, y );
	VecNorm( y, NORM_2, &v_norm );
	
	PetscPrintf( PETSC_COMM_WORLD, "norm(A {1}) = %5.5e \n", v_norm );
	PetscPrintf( PETSC_COMM_WORLD, "norm(A) = %5.5e \n\n", m_norm );
	
	
	/* Write matrix out in MM format */
	PetscPrintf( PETSC_COMM_WORLD, "Dump as MM \n");
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, "petsc.mtx", &viewer );
	__MatViewMatrixMarket( A, viewer );
	Stg_PetscViewerDestroy( & viewer );
	
	MatDestroy( & A );
	

	/* Read it back in and compute norms */
	__MatLoadMatrixMarket( PETSC_COMM_WORLD, "petsc.mtx", MATAIJ, &An );
	PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO );
	MatView( An, PETSC_VIEWER_STDOUT_WORLD );
	
	MatNorm( An, NORM_1, &m_norm );
	
	VecSet( y, 0.0 );
	VecSet( x, 1.0 );
	
	MatMult( An,x, y );
	VecNorm( y, NORM_2, &v_norm );
	
	PetscPrintf( PETSC_COMM_WORLD, "norm(An {1}) = %5.5e \n", v_norm );
	PetscPrintf( PETSC_COMM_WORLD, "norm(An) = %5.5e \n\n", m_norm );
	
	
	MatDestroy( & An );
	Stg_VecDestroy( & x );
	Stg_VecDestroy( & y );
}







int main( int argc, char **args )
{
	
	PetscInitialize( &argc, &args, (char *)0, PETSC_NULL );
	
	test_VecLoadMatrixMarket();
	test_MatAIJLoadMatrixMarket();
	
	PetscFinalize();
	return 0;
	
}
