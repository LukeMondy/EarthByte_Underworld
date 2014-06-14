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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>


#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "petscext_helpers.h"

PetscErrorCode VecLoad_MatrixMarket( MPI_Comm comm, const char fname[], VecType outtype, Vec *x )
{
	FILE *fp;
	PetscInt m,n, i, s,e,M, add;
	PetscScalar *_x;
	PetscInt *row_idx;
	size_t length;
	char ext[10];
	char header[2000];
	PetscTruth flg;
	//PetscInt ten_percent;
	PetscReal tmp;
	
	
	PetscFunctionBegin;
	PetscPrintf( comm, "VecCreateFromMatrixMarket: Reading file %s \n", fname );
	
	/* check extension */
	PetscStrlen( fname, &length );
	PetscStrcpy( ext, &fname[length-3] );
	PetscStrcmp( ext, "mtx", &flg );
	if( flg == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Extension must be .mtx" );
	}
	
	/* open file to read */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	
	/* read header */
	fgets( header, 2000, fp );
	PetscPrintf( comm, "  %s", header );
	/* read vector size, row x 1 */
	fscanf( fp, "%d %d", &m, &n );
	PetscPrintf( comm, "  m=%d : n=%d \n", m,n );
	
	
	VecCreate( comm, x );
	VecSetSizes( *x, PETSC_DECIDE, m );
	if( outtype != PETSC_NULL) { VecSetType( *x, outtype ); }
	VecSetFromOptions( *x );
	
	VecGetOwnershipRange( *x, &s, &e );
	VecGetLocalSize( *x, &M );
	
	
	/* allocate array of length m, read file contents into array */
	_x = (PetscScalar*)malloc( sizeof(PetscScalar) * M );
	row_idx = (PetscInt*)malloc( sizeof(PetscInt) * M );
	for( i=0; i<M; i++ ) {
		row_idx[i] = s + i;
	}
	
	
	/* Now only save the parts needed for the local part of vector */
	/* Could do some funky file skip but its not really necessary */
	//ten_percent =(PetscInt)(  10 * ((double)m / 100.0) );
	add = 0;
	for( i=0; i<m; i++ ) {
#ifdef PETSC_USE_COMPLEX
		PetscReal re,im;
#endif
		
		fscanf( fp, "%lf", &tmp );
		if( i >= s && i < e ) {
			_x[add] = tmp;
			
#ifdef PETSC_USE_COMPLEX
			re = tmp;
			im = 0.0;
			PetscRealPart(_x[add])      = re;
			PetscImaginaryPart(_x[add]) = im;
#else
			_x[add] = tmp;
#endif 
			
			add++;
		}
		/*
		if( i%ten_percent == 0 ) {
			PetscPrintf( comm, "\tInserted %d of %d entries \n", i, m );
		}
		*/
	}
	
	VecSetValues( *x, M, row_idx, _x, INSERT_VALUES );
	VecAssemblyBegin( *x );
	VecAssemblyEnd( *x );
	
	free( _x );
	free( row_idx );
	fclose( fp );
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecView_MatrixMarket( Vec x, PetscViewer viewer )
{
	int size;
	MPI_Comm comm;
	PetscInt N;
	
	
	PetscFunctionBegin;
	
	PetscObjectGetComm( (PetscObject)x, &comm );
	MPI_Comm_size( comm, &size );
	
	if( size != 1 ) {
		Vec seq_x;
		VecScatter ctx;
		
		VecGetSize( x, &N );
		VecCreate( PETSC_COMM_SELF, &seq_x );
		VecSetSizes( seq_x, PETSC_DECIDE, N );
		VecSetFromOptions( seq_x );
		
		VecScatterCreateToZero(x,&ctx,&seq_x);
		VecScatterBegin(ctx,x,seq_x,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(ctx,x,seq_x,INSERT_VALUES,SCATTER_FORWARD);
		Stg_VecScatterDestroy(&ctx);
		
		VecView_MatrixMarket( seq_x, viewer );
		Stg_VecDestroy( &seq_x );
	}
	else {
		PetscScalar *xa;
		PetscInt i;
		
		VecGetSize( x, &N );
		VecGetArray( x, &xa );
		
		PetscViewerASCIIPrintf( viewer, "%%%%MatrixMarket matrix array real general\n" );
		PetscViewerASCIIPrintf( viewer, "%d 1\n", N );
		for( i=0; i<N; i++ ) {
			PetscViewerASCIIPrintf( viewer, " %1.13e\n", xa[i] );
		}
		VecRestoreArray( x, &xa );
		
	}
	PetscFunctionReturn(0);
}




/* DEPRECIATED */
/*
Scan file twice.
First time we compute the number of nonzeros per row.
Second time we assemble the matrix.
*/
PetscErrorCode MatSEQAIJCreateFromMatrixMarket( MPI_Comm comm, const char fname[], Mat *A )
{
	PetscErrorCode ierr;
	FILE *fp;
	PetscInt m,n,nnz, k,i,j;
	PetscReal val;
	size_t length;
	char ext[10];
	char header[2000];
	PetscInt *nz;//, ten_percent;
	PetscTruth flg;
#ifdef PETSC_USE_COMPLEX
	PetscReal re, im;
#endif
	PetscScalar scalar;
	
	
	PetscFunctionBegin;
	PetscPrintf( comm, "MatSEQAIJCreateFromMatrixMarket: Reading file %s \n", fname );
	
	/* check extension */
	PetscStrlen( fname, &length );
	PetscStrcpy( ext, &fname[length-3] );
	PetscStrcmp( ext, "mtx", &flg );
	if( flg == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Extension must be .mtx" );
	}
	
	/* PHASE 1. Determine nonzero count per row */
	/* open file to read */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	/* read header */
	fgets( header, 2000, fp );
	PetscPrintf( comm, "  %s", header );
	
	/* read matrix size, row col row*col */
	fscanf( fp, "%d %d %d", &m, &n, &nnz );
	//mn = m * n;
	PetscPrintf( comm, "  m=%d : n=%d : nnz=%d\n", m,n,nnz );
	
	
	MatCreate( comm, A );
	MatSetSizes( *A, PETSC_DECIDE,PETSC_DECIDE, m,n );
	MatSetType( *A, MATSEQAIJ );
	
	
	/* Determine nonzero structure */
	nz = (PetscInt*)malloc( sizeof(PetscInt)*m );
	for( i=0; i<m; i++ ) nz[i] = 0;
	
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
		nz[ i-1 ]++;
	}
	fclose( fp );
	
	MatSeqAIJSetPreallocation( *A, PETSC_NULL, nz );
	
	
	/* PHASE 2. Assemble matrix */
	/* open file to read again */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	/* read header */
	fgets( header, 2000, fp );
	fscanf( fp, "%d %d %d", &m, &n, &nnz );
	
	//ten_percent = (PetscInt)( 10 * ((double)nnz / 100.0) );
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
	//	printf("%d %d %14.13e \n", i,j,val );
		
#ifdef PETSC_USE_COMPLEX
		re = val;
		im = 0.0;
		PetscRealPart(scalar)      = re;
		PetscImaginaryPart(scalar) = im;
#else
		scalar = val;
#endif 
		
		
		ierr=MatSetValue( *A, i-1,j-1, scalar, INSERT_VALUES );CHKERRQ(ierr);
		/*
		if( k%ten_percent == 0 ) {
			PetscPrintf( comm, "\tInserted %d of %d entries\n", k, nnz );
		}
		*/
	}
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
	
	
	free( nz );
	fclose(fp);
	
	PetscFunctionReturn(0);
}


PetscErrorCode MatMPIAIJCreateFromMatrixMarket( MPI_Comm comm, const char fname[], Mat *A )
{
	PetscErrorCode ierr;
	FILE *fp;
	PetscInt m,n, M,N, nnz, k,i,j, s,e;
	PetscReal val;
	size_t length;
	char ext[10];
	char header[2000];
	PetscInt *dnz, *onz, row_index, col_index;
	//PetscInt ten_percent;
	PetscTruth flg;
#ifdef PETSC_USE_COMPLEX
	PetscReal re,im;
#endif
	PetscScalar scalar;
	
	
	PetscFunctionBegin;
	PetscPrintf( comm, "MatMPIAIJCreateFromMatrixMarket: Reading file %s \n", fname );
	
	/* check extension */
	PetscStrlen( fname, &length );
	PetscStrcpy( ext, &fname[length-3] );
	PetscStrcmp( ext, "mtx", &flg );
	if( flg == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Extension must be .mtx" );
	}
	
	/* PHASE 1. Determine nonzero count per row */
	/* open file to read */
	fp = fopen( fname, "r" );
	if (!fp) {	Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );	}
	
	/* read header */
	fgets( header, 2000, fp );
	PetscPrintf( comm, "  %s", header );
	
	/* read matrix size, row col row*col */
	fscanf( fp, "%d %d %d", &M, &N, &nnz );
	//mn = M * N;
	PetscPrintf( comm, "  m=%d : n=%d : nnz=%d\n", M,N,nnz );
	
	
	MatCreate( comm, A );
	MatSetSizes( *A, PETSC_DECIDE,PETSC_DECIDE, M,N );
	MatSetType( *A, MATAIJ );
	
	/* Calling ownership range will compute the local m,n if there were given by PETSC_DECIDE in MatSetSizes() */
	MatGetOwnershipRange( *A, &s, &e );
	MatGetLocalSize( *A, &m, &n );
	
	/* Determine nonzero structure */
	dnz = (PetscInt*)malloc( sizeof(PetscInt)*m );
	onz = (PetscInt*)malloc( sizeof(PetscInt)*m );
	
	/* init */
	for( i=0; i<m; i++ ) {
		dnz[i] = 0;
		onz[i] = 0;
	}
	
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
		row_index = i-1;
		col_index = j-1;
		
		if( row_index >= s && row_index < e ) {
			if( col_index >= s && col_index < e ) {
				dnz[ row_index - s ]++;    /* shift by -s to make global row index a local row index */
			}
			else {
				onz[ row_index - s ]++;    /* shift by -s to make global row index a local row index */
			}
		}
	}
	fclose( fp );
	
	MatMPIAIJSetPreallocation( *A, PETSC_NULL,dnz, PETSC_NULL,onz );
	
	
	/* PHASE 2. Assemble matrix */
	/* open file to read again */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	/* read header */
	fgets( header, 2000, fp );
	fscanf( fp, "%d %d %d", &M, &N, &nnz );
	//ten_percent = (PetscInt)( 10 * ((double)nnz / 100.0) );
	
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
		row_index = i - 1;
		col_index = j - 1;
		
#ifdef PETSC_USE_COMPLEX
		re = val;
		im = 0.0;
		PetscRealPart(scalar)      = re;
		PetscImaginaryPart(scalar) = im;
		
		if( row_index >= s && row_index < e ) {
			ierr=MatSetValue( *A, row_index, col_index, scalar, INSERT_VALUES );CHKERRQ(ierr);
		}
#else
		scalar = val;
		
		if( row_index >= s && row_index < e ) {
			ierr=MatSetValue( *A, row_index, col_index, scalar, INSERT_VALUES );CHKERRQ(ierr);
		}
#endif 
		
		
		
		/*
		if( k%ten_percent == 0 ) {
			PetscPrintf( comm, "\tInserted %d of %d entries\n", k, nnz );
		}
		*/
	}
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
	
	
	free( dnz );
	free( onz );
	fclose(fp);
	
	PetscFunctionReturn(0);
}


/*
Scan file twice.
First time we compute the number of nonzeros per row.
Second time we assemble the matrix.
*/
PetscErrorCode MatSeqAIJCreateBinaryFromMatrixMarket( const char fname[], const char fname_out[] )
{
	PetscErrorCode ierr;
	FILE *fp;
	PetscInt m,n, nnz, k,i,j;
	PetscReal val;
	size_t length;
	char ext[10];
	char header[2000];
	PetscInt *nz;//, ten_percent;
	PetscTruth flg;
#ifdef PETSC_USE_COMPLEX
	PetscReal re, im;
#endif
	PetscScalar scalar;
	Mat seq_A;
	PetscViewer viewer;
	
	
	
	PetscFunctionBegin;
	PetscPrintf( PETSC_COMM_SELF, "MatSeqAIJCreateBinaryFromMatrixMarket: Reading file %s \n", fname );
	
	
	/* check extension */
	PetscStrlen( fname, &length );
	PetscStrcpy( ext, &fname[length-3] );
	PetscStrcmp( ext, "mtx", &flg );
	if( flg == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Extension must be .mtx" );
	}
	
	/* PHASE 1. Determine nonzero count per row */
	/* open file to read */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	/* read header */
	fgets( header, 2000, fp );
	PetscPrintf( PETSC_COMM_SELF, "  %s", header );
	
	/* read matrix size, row col row*col */
	fscanf( fp, "%d %d %d", &m, &n, &nnz );
	//mn = m * n;
	PetscPrintf( PETSC_COMM_SELF, "  m=%d : n=%d : nnz=%d\n", m,n,nnz );
	
	
	MatCreate( PETSC_COMM_SELF, &seq_A );
	MatSetSizes( seq_A, PETSC_DECIDE,PETSC_DECIDE, m,n );
	MatSetType( seq_A, MATSEQAIJ );
	
	
	/* Determine nonzero structure */
	nz = (PetscInt*)malloc( sizeof(PetscInt)*m );
	for( i=0; i<m; i++ ) nz[i] = 0;
	
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
		nz[ i-1 ]++;
	}
	fclose( fp );
	
	MatSeqAIJSetPreallocation( seq_A, PETSC_NULL, nz );
	
	
	/* PHASE 2. Assemble matrix */
	/* open file to read again */
	fp = fopen( fname, "r" );
	if( fp == PETSC_NULL ) {
		Stg_SETERRQ1( PETSC_ERR_SUP, "Cannot open file %s", fname );
	}
	
	/* read header */
	fgets( header, 2000, fp );
	fscanf( fp, "%d %d %d", &m, &n, &nnz );
	
	//ten_percent = (PetscInt)( 10 * ((double)nnz / 100.0) );
	for( k=0; k<nnz; k++ ) {
		fscanf( fp, "%d %d %lf", &i, &j, &val );
	//	printf("%d %d %14.13e \n", i,j,val );
		
#ifdef PETSC_USE_COMPLEX
		re = val;
		im = 0.0;
		PetscRealPart(scalar)      = re;
		PetscImaginaryPart(scalar) = im;
#else
		scalar = val;
#endif 
		
		
		ierr=MatSetValue( seq_A, i-1,j-1, scalar, INSERT_VALUES );CHKERRQ(ierr);
		/*
		if( k%ten_percent == 0 ) {
			PetscPrintf( PETSC_COMM_SELF, "\tInserted %d of %d entries\n", k, nnz );
		}
		*/
	}
	MatAssemblyBegin(seq_A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(seq_A,MAT_FINAL_ASSEMBLY);
	
	
	free( nz );
	fclose(fp);
	
	
	/* Now we have fully assembled sequential matrix */
	/* Create viewer and write to binary */
//	PetscViewerCreate(PETSC_COMM_SELF,&viewer);
//	PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE);
//	PetscViewerFileSetName(viewer, fname_out);
	PetscViewerBinaryOpen(PETSC_COMM_SELF,fname_out,FILE_MODE_WRITE,&viewer);
	MatView(seq_A,viewer);
	
	Stg_PetscViewerDestroy(&viewer);
	Stg_MatDestroy(&seq_A);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode __Stg_MatLoadMatrixMarket( MPI_Comm comm, const char fname[], const MatType type, Mat *A )
{
	PetscViewer viewer;
	time_t      currTime;
	struct tm*  timeInfo;
	int         adjustedYear;
	int         adjustedMonth;
	char        *fname_out;
	int z_adjustedYear, z_adjustedMonth, z_tm_mday, z_tm_hour, z_tm_min, z_tm_sec;
	PetscMPIInt rank;
	
	PetscFunctionBegin;
	
	
	/* get file name with date */
	currTime = time( NULL );
	timeInfo = localtime( &currTime );
	/* See man localtime() for why to adjust these */
	adjustedYear = 1900 + timeInfo->tm_year;
	adjustedMonth = 1 + timeInfo->tm_mon;

	z_adjustedYear = adjustedYear;
	z_adjustedMonth = adjustedMonth;
	z_tm_mday = timeInfo->tm_mday;
	z_tm_hour = timeInfo->tm_hour;
	z_tm_min = timeInfo->tm_min;
	z_tm_sec = timeInfo->tm_sec;
	
	MPI_Bcast ( &z_adjustedYear, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_adjustedMonth, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_mday, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_hour, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_min, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_sec, 1, MPI_INT, 0, comm );

	/* Format; tmp-A-Stg_MatLoadMM--YYYY.MM.DD-HH.MM.SS.mtx */	
	asprintf( &fname_out, "tmp-A-Stg_MatLoadMM--%.4d.%.2d.%.2d-%.2d.%.2d.%.2d.mtx",
			z_adjustedYear, z_adjustedMonth, z_tm_mday,
			z_tm_hour, z_tm_min, z_tm_sec );
	
	MPI_Comm_rank( comm, &rank );
	if(rank==0) {
		MatSeqAIJCreateBinaryFromMatrixMarket( fname, fname_out );
	}
	/* make sure proc zero has finished writing before proceeding */
	MPI_Barrier(comm);
	
	PetscViewerBinaryOpen(comm,fname_out,FILE_MODE_READ,&viewer);
	Stg_MatLoad(viewer, type, A);
	
	Stg_PetscViewerDestroy(&viewer);
	free(fname_out);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode MatAIJLoad_MatrixMarket( MPI_Comm comm, const char fname[], Mat *A )
{
	int size;
	
	PetscFunctionBegin;
	
	MPI_Comm_size( comm, &size );
	if( size != 1 ) {
		MatMPIAIJCreateFromMatrixMarket( comm, fname, A );
	}
	else {
		MatSEQAIJCreateFromMatrixMarket( comm, fname, A );
	}
	
	PetscFunctionReturn(0);
}


PetscErrorCode MatView_MatrixMarket_seq( Mat A, PetscViewer viewer )
{
	PetscInt ridx,i,j;
	PetscInt nc, cidx, M,N, m,n, s,e;
	const PetscInt *colidx;
	const PetscScalar *vals;
	MatInfo mat_info;
	PetscLogDouble NNZ;
	
	
	PetscFunctionBegin;
	
	MatGetSize( A, &M,&N );
	MatGetLocalSize( A, &m,&n );
	MatGetOwnershipRange( A, &s, &e );
	
	/*
	%%MatrixMarket matrix coordinate real general
	M N NNZ
	*/
	
	// MAT_GLOBAL_MAX or MAT_GLOBAL_SUM
	MatGetInfo( A, MAT_GLOBAL_SUM, &mat_info );
	NNZ = mat_info.nz_used;
	
	PetscViewerASCIIPrintf( viewer, "%%%%MatrixMarket matrix coordinate real general\n" );
	PetscViewerASCIIPrintf( viewer, "%d %d %d\n", M,N, (PetscInt)NNZ );
	
	for( i=0; i<m; i++ ) {
		ridx = i+s;
		MatGetRow( A, ridx, &nc, &colidx, &vals );
		
		for( j=0; j<nc; j++ ) {
			cidx = colidx[j];
			PetscViewerASCIIPrintf( viewer, "%d %d %1.13e\n", ridx+1,cidx+1, vals[j] );
		}
		
		MatRestoreRow( A, ridx, &nc, &colidx, &vals );
	}
	
	PetscFunctionReturn(0);
}

/*
This is pretty crappy due to the massive amount of redundnent work being done.
We create a matrix with identical entries on each processor. Proc 0 dumps the info.
We generate these redundent matrices once for each processor.
*/
PetscErrorCode MatView_MatrixMarket_mpi( Mat A, PetscViewer viewer )
{
	PetscInt ridx,i,j;
	PetscInt nc, cidx, M,N, m,n, s,e;
	const PetscInt *colidx;
	const PetscScalar *vals;
	MatInfo mat_info;
	PetscLogDouble NNZ;
	int size,rank,p;
	MPI_Comm comm;
	PetscInt *ranges;
	IS isrow, iscol;
	Mat newmat;
	
	
	PetscFunctionBegin;
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MPI_Comm_size( comm, &size );
	MPI_Comm_rank( comm, &rank );
	
	MatGetOwnershipRanges( A ,(const PetscInt**)&ranges );
	
	MatGetSize( A, &M,&N );
	MatGetLocalSize( A, &m,&n );
	MatGetOwnershipRange( A, &s, &e );
	
	/*
	%%MatrixMarket matrix coordinate real general
	M N NNZ
	*/
	
	// MAT_GLOBAL_MAX or MAT_GLOBAL_SUM
	MatGetInfo( A, MAT_GLOBAL_SUM, &mat_info );
	NNZ = mat_info.nz_used;
	
	PetscViewerASCIIPrintf( viewer, "%%%%MatrixMarket matrix coordinate real general\n" );
	PetscViewerASCIIPrintf( viewer, "%d %d %d\n", M,N, (PetscInt)NNZ );
	
	
	ISCreateStride( comm, N,0,1, &iscol );
	for( p=0; p<size; p++ ) {
		
		ISCreateStride( comm, (ranges[p+1]-ranges[p]),ranges[p],1, &isrow );
		MatGetSubMatrix( A,isrow,iscol, MAT_INITIAL_MATRIX, &newmat );
		
		if( rank == 0 ) {
			for( i=0; i<(ranges[p+1]-ranges[p]); i++ ) {
				
				ridx = i+ranges[p];
				MatGetRow( newmat, i, &nc, &colidx, &vals );
				
				for( j=0; j<nc; j++ ) {
					cidx = colidx[j];
					PetscViewerASCIIPrintf( viewer, "%d %d %1.13e\n", ridx+1,cidx+1, vals[j] );
				}
				
				MatRestoreRow( newmat, i, &nc, &colidx, &vals );
			}
		}
		
		
		Stg_MatDestroy( &newmat );
		Stg_ISDestroy( &isrow );
	}
	Stg_ISDestroy( &iscol );
	
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode MatView_MatrixMarket( Mat A, PetscViewer viewer )
{
	int size;
	MPI_Comm comm;
	
	
	PetscFunctionBegin;
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MPI_Comm_size( comm, &size );
	
	if( size != 1 ) {
		MatView_MatrixMarket_mpi( A, viewer );
	}
	else {
		MatView_MatrixMarket_seq( A, viewer );
	}
	PetscFunctionReturn(0);
}



PetscErrorCode MatViewMatrixMarket_SEQAIJ( Mat A, PetscViewer viewer )
{
	PetscInt ridx,i,j;
	PetscInt nc, cidx, M,N, m,n, s,e;
	const PetscInt *colidx;
	const PetscScalar *vals;
	MatInfo mat_info;
	PetscLogDouble NNZ;
	
	
	PetscFunctionBegin;
	
	MatGetSize( A, &M,&N );
	MatGetLocalSize( A, &m,&n );
	MatGetOwnershipRange( A, &s, &e );
	
	/*
	%%MatrixMarket matrix coordinate real general
	M N NNZ
	*/
	
	// MAT_GLOBAL_MAX or MAT_GLOBAL_SUM
	MatGetInfo( A, MAT_GLOBAL_SUM, &mat_info );
	NNZ = mat_info.nz_used;
	
	PetscViewerASCIIPrintf( viewer, "%%%%MatrixMarket matrix coordinate real general\n" );
	PetscViewerASCIIPrintf( viewer, "%d %d %d\n", M,N, (PetscInt)NNZ );
	
	for( i=0; i<m; i++ ) {
		ridx = i+s;
		MatGetRow( A, ridx, &nc, &colidx, &vals );
		
		for( j=0; j<nc; j++ ) {
			cidx = colidx[j];
			PetscViewerASCIIPrintf( viewer, "%d %d %1.13e\n", ridx+1,cidx+1, vals[j] );
		}
		
		MatRestoreRow( A, ridx, &nc, &colidx, &vals );
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode __MatViewMatrixMarket( Mat A, PetscViewer mm_viewer )
{
	PetscMPIInt rank;
	MPI_Comm    comm;
	PetscViewer binary_viewer;
	PetscTruth  is_ascii;
	time_t      currTime;
	struct tm*  timeInfo;
	int         adjustedYear;
	int         adjustedMonth;
	char        *fname_out;
	int z_adjustedYear, z_adjustedMonth, z_tm_mday, z_tm_hour, z_tm_min, z_tm_sec;
	Mat seq_A;
	
	
	PetscFunctionBegin;
	
	/* check mm_viewer is either PETSC_VIEWER_STDOUT_WORLD or PETSC_VIEWER_STDOUT_SELF or PETSC_VIEWER_ASCII_MATLAB */
	PetscTypeCompare((PetscObject)mm_viewer,PETSC_VIEWER_ASCII,&is_ascii);
	if(is_ascii==PETSC_FALSE) {
		Stg_SETERRQ(PETSC_ERR_SUP,"Only for viewers of type PETSC_VIEWER_ASCII_MATLAB");
	}
	
	PetscObjectGetComm( (PetscObject)A, &comm );
	MPI_Comm_rank( comm, &rank );
	
	
	/* get file name with date */
	currTime = time( NULL );
	timeInfo = localtime( &currTime );
	/* See man localtime() for why to adjust these */
	adjustedYear = 1900 + timeInfo->tm_year;
	adjustedMonth = 1 + timeInfo->tm_mon;
	
	z_adjustedYear = adjustedYear;
	z_adjustedMonth = adjustedMonth;
	z_tm_mday = timeInfo->tm_mday;
	z_tm_hour = timeInfo->tm_hour;
	z_tm_min = timeInfo->tm_min;
	z_tm_sec = timeInfo->tm_sec;
	
	MPI_Bcast ( &z_adjustedYear, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_adjustedMonth, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_mday, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_hour, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_min, 1, MPI_INT, 0, comm );
	MPI_Bcast ( &z_tm_sec, 1, MPI_INT, 0, comm );
	
	/* Format; tmp-A-MatViewMM--YYYY.MM.DD-HH.MM.SS.mtx */
	asprintf( &fname_out, "tmp-A-MatViewMM--%.4d.%.2d.%.2d-%.2d.%.2d.%.2d.mtx",
			z_adjustedYear, z_adjustedMonth, z_tm_mday,
			z_tm_hour, z_tm_min, z_tm_sec );

	
	PetscViewerBinaryOpen(comm,fname_out,FILE_MODE_WRITE,&binary_viewer);
	MatView(A,binary_viewer);
	Stg_PetscViewerDestroy(&binary_viewer);
	/* make sure all procs have finished writing before proceeding */
	MPI_Barrier(comm);
	
	PetscViewerBinaryOpen(PETSC_COMM_SELF,fname_out,FILE_MODE_READ,&binary_viewer);
	Stg_MatLoad(binary_viewer,MATSEQAIJ,&seq_A);
	if( rank == 0 ) {
		MatViewMatrixMarket_SEQAIJ( seq_A, mm_viewer );
	}
	Stg_MatDestroy(&seq_A);
	Stg_PetscViewerDestroy(&binary_viewer);
	
	free(fname_out);
	
	PetscFunctionReturn(0);
}

