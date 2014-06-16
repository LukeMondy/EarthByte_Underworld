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


static char help[] = "Utility to convert petsc matrices(vectors) from the native binaary representation to an ascii MatrixMarket format. \n\
Input arguments are: \n\
    -mat /path/to/petsc.mat  (to convert a matrix)\n\
    -vec /path/to/petsc.vec  (to convert a vector)\n\
    -output output_filename  (to specifiy the name of the converted object, defaults used are matrix.mtx or vector.mtx) \n\n";


void convert_vec( const char infile[], const char outfile[] )
{
	Vec x;
	PetscViewer viewer;
	
	/* Read petsc binary */
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, infile, FILE_MODE_READ, &viewer );
	Stg_VecLoad( viewer, PETSC_NULL, &x );
	Stg_PetscViewerDestroy( &viewer );
	
	/* Write Matrix Market */
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, outfile, &viewer );
	VecView_MatrixMarket( x, viewer );
	Stg_PetscViewerDestroy( & viewer );
	
	Stg_VecDestroy( & x );
}

void convert_mat( const char infile[], const char outfile[] )
{
	Mat A;
	PetscViewer viewer;
	
	/* Read petsc binary */
	PetscViewerBinaryOpen( PETSC_COMM_WORLD, infile, FILE_MODE_READ, &viewer );
	Stg_MatLoad( viewer, MATSEQAIJ, &A );
	Stg_PetscViewerDestroy( & viewer );
	
	/* Write Matrix Market */
	PetscViewerASCIIOpen( PETSC_COMM_WORLD, outfile, &viewer );
	MatView_MatrixMarket( A, viewer );
	Stg_PetscViewerDestroy( & viewer );
	
	MatDestroy( & A );
}


int main( int argc, char **args )
{
	PetscMPIInt size;
	PetscTruth mat_flg, vec_flg, out_flg;
	char infile[PETSC_MAX_PATH_LEN];
	char outfile[PETSC_MAX_PATH_LEN];
	size_t len;
	
	
	PetscInitialize( &argc, &args, (char *)0, help );
	MPI_Comm_size( PETSC_COMM_WORLD, &size );
	if( size != 1 ) Stg_SETERRQ( PETSC_ERR_SUP, "Uniprocessor only" );
	
	len = PETSC_MAX_PATH_LEN -1;
	mat_flg = vec_flg = PETSC_FALSE;
	PetscOptionsGetString( PETSC_NULL, "-mat", infile, len, &mat_flg );
	PetscOptionsGetString( PETSC_NULL, "-vec", infile, len, &vec_flg );
	
	PetscOptionsGetString( PETSC_NULL, "-output", outfile, len, &out_flg );
	if( !out_flg ) {
		if( mat_flg ) {
			sprintf( outfile, "matrix.mtx" );
		}
		if( vec_flg ) {
			sprintf( outfile, "vec.mtx" );
		}
	}
	
	if( mat_flg ) {
		PetscPrintf( PETSC_COMM_WORLD, "Input matrix: %s \n", infile );
		PetscPrintf( PETSC_COMM_WORLD, "Output matrix: %s \n", outfile );
		
		convert_mat( infile, outfile );
	}
	if( vec_flg ) {
		PetscPrintf( PETSC_COMM_WORLD, "Input vector: %s \n", infile );
		PetscPrintf( PETSC_COMM_WORLD, "Output vector: %s \n", outfile );
		
		convert_vec( infile, outfile );
	}
	if( !mat_flg && !vec_flg ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "You must specify a matrix (or vector) to convert via -mat xxx (-vec xxx)" );
	}
	
	PetscFinalize();
	return 0;
	
}
