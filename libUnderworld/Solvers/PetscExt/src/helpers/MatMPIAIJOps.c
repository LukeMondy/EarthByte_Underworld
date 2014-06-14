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

#include "private/matimpl.h"
#include "src/mat/impls/aij/mpi/mpiaij.h"

#include "private/compat/petsccompat.h"

PetscErrorCode PetscExtMatGetRowMax_MPIAIJ(Mat mat,Vec v,PetscInt idx[])
{
	Mat_MPIAIJ     *aij = (Mat_MPIAIJ*)mat->data;
	Mat            A = aij->A;
	Mat            B = aij->B;
	PetscInt       i,mA,nA, mB,nB;
	Vec            vA,vB;
	PetscScalar    *_vA, *_vB, *_v;
	PetscInt       *idxA,*idxB;
	PetscScalar    a,b;
	MatInfo        info;
	
	MatGetLocalSize( A, &mA, &nA );
	MatGetLocalSize( B, &mB, &nB );
	
	MatGetVecs( A, PETSC_NULL, &vA );
	MatGetVecs( B, PETSC_NULL, &vB );
	VecSet( vA, 0.0 );
	VecSet( vB, 0.0 );
	
	if (idx!=PETSC_NULL) { 
		PetscMalloc( sizeof(PetscInt) * mA, &idxA ); 
		PetscMalloc( sizeof(PetscInt) * mB, &idxB ); 
	}
	else {
		idxA = idxB = PETSC_NULL;
	}
	
	MatGetInfo( A, MAT_GLOBAL_MAX, &info );
	if( info.nz_allocated > 0.0 ) {   MatGetRowMax( A, vA, idxA );   }

	MatGetInfo( B, MAT_GLOBAL_MAX, &info );
	if( info.nz_allocated > 0.0 ) {   MatGetRowMax( B, vB, idxB );   }
	
	VecGetArray( vA, &_vA );
	VecGetArray( vB, &_vB );
	VecGetArray( v, &_v );
	
	for (i=0; i<mA; i++) {
		a = _vA[i];
		b = _vB[i];
		
		/* compare  (a > b) */
		if( (PetscRealPart(a)) > (PetscRealPart(b)) ) { //
			_v[i] = _vA[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxA[i];		}
		} //
		else { 
			_v[i] = _vB[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxB[i];		}
		}
	}
	
	
	/* tidy up */
	if (idx!=PETSC_NULL) {
		PetscFree( idxA );
		PetscFree( idxB );
	}
	
	VecRestoreArray( vA, &_vA );
	VecRestoreArray( vB, &_vB );
	VecRestoreArray( v, &_v );
	
	Stg_VecDestroy( &vA );
	Stg_VecDestroy( &vB );
	
	PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMax_MATAIJ(Mat A,Vec r,PetscInt idx[])
{
        PetscInt i,j;
        PetscInt m, n;

        const PetscInt *cols;
        PetscInt ncols,r_col;
        const PetscScalar *vals;
        PetscScalar *r_vals;
        PetscReal cmp, abv;


        MatGetOwnershipRange( A, &m, &n );

        VecGetArray( r, &r_vals );
        for( i=m; i<n; i++ ) {
                MatGetRow( A, i, &ncols, &cols, &vals );
		if( ncols==0 ) continue;
		
                cmp = PetscRealPart(vals[0]);
		r_col = cols[0];
                for( j=1; j<ncols; j++ ) {
			abv = PetscRealPart(vals[j]);
                        if( abv > cmp ) {  /* Look for maximum */
                                cmp = abv;
                                r_col = cols[j];
                        }
                }
		if( idx ) {  idx[i-m] = r_col;  }
                r_vals[i-m] = cmp;

                MatRestoreRow( A, i, &ncols, &cols, &vals );
        }
        VecRestoreArray( r, &r_vals );

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMax(Mat mat,Vec vec,PetscInt idx[])
{
        PetscTruth is_mpiaij;

        PetscTypeCompare( (PetscObject)mat, MATMPIAIJ, &is_mpiaij );
        if (is_mpiaij==PETSC_TRUE) {    PetscExtMatGetRowMax_MATAIJ(mat,vec,idx);       }
        else {                          MatGetRowMax(mat,vec,idx);                      }

        PetscFunctionReturn(0);
}


PetscErrorCode PetscExtMatGetRowMin_MPIAIJ(Mat mat,Vec v,PetscInt idx[])
{
	Mat_MPIAIJ     *aij = (Mat_MPIAIJ*)mat->data;
	Mat            A = aij->A;
	Mat            B = aij->B;
	PetscInt       i,mA,nA, mB,nB;
	Vec            vA,vB;
	PetscScalar    *_vA, *_vB, *_v;
	PetscInt       *idxA,*idxB;
	PetscScalar    a,b;
	MatInfo        info;	

	
	MatGetLocalSize( A, &mA, &nA );
	MatGetLocalSize( B, &mB, &nB );
	
	MatGetVecs( A, PETSC_NULL, &vA );
	MatGetVecs( B, PETSC_NULL, &vB );
	VecSet( vA, 0.0 );
	VecSet( vB, 0.0 );

	if (idx!=PETSC_NULL) { 
		PetscMalloc( sizeof(PetscInt) * mA, &idxA ); 
		PetscMalloc( sizeof(PetscInt) * mB, &idxB ); 
	}
	else {
		idxA = idxB = PETSC_NULL;
	}
	
        MatGetInfo( A, MAT_GLOBAL_MAX, &info );
        if( info.nz_allocated > 0.0 ) {   MatGetRowMin( A, vA, idxA );   }

        MatGetInfo( B, MAT_GLOBAL_MAX, &info );
        if( info.nz_allocated > 0.0 ) {   MatGetRowMin( B, vB, idxB );   }
	
	VecGetArray( vA, &_vA );
	VecGetArray( vB, &_vB );
	VecGetArray( v, &_v );
	
	for (i=0; i<mA; i++) {
		a = _vA[i];
		b = _vB[i];
		
		/* compare  (a < b) */
		if( (PetscRealPart(a)) < (PetscRealPart(b)) ) { //
			_v[i] = _vA[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxA[i];		}
		} //
		else { 
			_v[i] = _vB[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxB[i];		}
		}
	}
	
	
	/* tidy up */
	if (idx!=PETSC_NULL) {
		PetscFree( idxA );
		PetscFree( idxB );
	}
	
	VecRestoreArray( vA, &_vA );
	VecRestoreArray( vB, &_vB );
	VecRestoreArray( v, &_v );
	
	Stg_VecDestroy( &vA );
	Stg_VecDestroy( &vB );
	
	PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMin_MATAIJ(Mat A,Vec r,PetscInt idx[])
{
        PetscInt i,j;
        PetscInt m, n;

        const PetscInt *cols;
        PetscInt ncols,r_col;
        const PetscScalar *vals;
        PetscScalar *r_vals;
        PetscReal cmp, abv;


        MatGetOwnershipRange( A, &m, &n );

        VecGetArray( r, &r_vals );
        for( i=m; i<n; i++ ) {
                MatGetRow( A, i, &ncols, &cols, &vals );
                if( ncols==0 ) continue;

                cmp = PetscRealPart(vals[0]);
                r_col = cols[0];
                for( j=1; j<ncols; j++ ) {
                        abv = PetscRealPart(vals[j]);
                        if( abv < cmp ) {  /* Look for minimum */
                                cmp = abv;
                                r_col = cols[j];
                        }
                }
		if( idx ) {  idx[i-m] = r_col; }
                r_vals[i-m] = cmp;

                MatRestoreRow( A, i, &ncols, &cols, &vals );
        }
        VecRestoreArray( r, &r_vals );

        PetscFunctionReturn(0);
}


PetscErrorCode PetscExtMatGetRowMin(Mat mat,Vec vec,PetscInt idx[])
{       
        PetscTruth is_mpiaij;

        PetscTypeCompare( (PetscObject)mat, MATMPIAIJ, &is_mpiaij );
        if (is_mpiaij==PETSC_TRUE) {    PetscExtMatGetRowMin_MATAIJ(mat,vec,idx);       }
        else {                          MatGetRowMin(mat,vec,idx);                      }

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMaxAbs_MPIAIJ(Mat mat,Vec v,PetscInt idx[])
{
	Mat_MPIAIJ     *aij = (Mat_MPIAIJ*)mat->data;
	Mat            A = aij->A;
	Mat            B = aij->B;
	PetscInt       i,mA,nA, mB,nB;
	Vec            vA,vB;
	PetscScalar    *_vA, *_vB, *_v;
	PetscInt       *idxA,*idxB;
	PetscScalar    a,b;
	PetscReal      ta,tb;
	MatInfo        info;	
	
	MatGetLocalSize( A, &mA, &nA );
	MatGetLocalSize( B, &mB, &nB );
	
	MatGetVecs( A, PETSC_NULL, &vA );
	MatGetVecs( B, PETSC_NULL, &vB );
	VecSet( vA, 0.0 );
	VecSet( vB, 0.0 );	


	if (idx!=PETSC_NULL) { 
		PetscMalloc( sizeof(PetscInt) * mA, &idxA ); 
		PetscMalloc( sizeof(PetscInt) * mB, &idxB ); 
	}
	else {
		idxA = idxB = PETSC_NULL;
	}

        MatGetInfo( A, MAT_GLOBAL_MAX, &info );	
        if( info.nz_allocated > 0.0 ) {   MatGetRowMaxAbs( A, vA, idxA );   }

        MatGetInfo( B, MAT_GLOBAL_MAX, &info );
        if( info.nz_allocated > 0.0 ) {   MatGetRowMaxAbs( B, vB, idxB );   }	


	VecGetArray( vA, &_vA );
	VecGetArray( vB, &_vB );
	VecGetArray( v, &_v );
	
	for (i=0; i<mA; i++) {
		a = _vA[i];
		b = _vB[i];
		
		ta = PetscAbsScalar(a);
		tb = PetscAbsScalar(b);
		
		/* compare  (ta > tb) */
		if (ta > tb) { //
			_v[i] = _vA[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxA[i];		}
		} //
		else { 
			_v[i] = _vB[i];
			if (idx!=PETSC_NULL) {		idx[i] = idxB[i];		}
		}
	}
	
	
	/* tidy up */
	if (idx!=PETSC_NULL) {
		PetscFree( idxA );
		PetscFree( idxB );
	}
	
	VecRestoreArray( vA, &_vA );
	VecRestoreArray( vB, &_vB );
	VecRestoreArray( v, &_v );
	
	Stg_VecDestroy( &vA );
	Stg_VecDestroy( &vB );
	
	PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMaxAbs_MATAIJ(Mat A,Vec r,PetscInt idx[])
{
        PetscInt i,j;
        PetscInt m, n;

        const PetscInt *cols;
        PetscInt ncols,r_col;
        const PetscScalar *vals;
        PetscScalar *r_vals;
        PetscReal cmp, abv;


        MatGetOwnershipRange( A, &m, &n );

        VecGetArray( r, &r_vals );
        for( i=m; i<n; i++ ) {
                MatGetRow( A, i, &ncols, &cols, &vals );
                if( ncols==0 ) continue;

                cmp = PetscAbsScalar(vals[0]);
                r_col = cols[0];
                for( j=1; j<ncols; j++ ) {
                        abv = PetscAbsScalar(vals[j]);
                        if( abv > cmp ) {  /* Look for maximum */
                                cmp = abv;
                                r_col = cols[j];  
                        }
                }
                if( idx ) {  idx[i-m] = r_col;  }
                r_vals[i-m] = cmp;

                MatRestoreRow( A, i, &ncols, &cols, &vals );
        }
        VecRestoreArray( r, &r_vals );

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMaxAbs(Mat mat,Vec vec,PetscInt idx[])
{
        PetscTruth is_mpiaij;

        PetscTypeCompare( (PetscObject)mat, MATMPIAIJ, &is_mpiaij );
        if (is_mpiaij==PETSC_TRUE) {    PetscExtMatGetRowMaxAbs_MATAIJ(mat,vec,idx);       }
        else {                          MatGetRowMaxAbs(mat,vec,idx);                      }

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMinAbs_MPIAIJ(Mat mat,Vec v,PetscInt idx[])
{
        Mat_MPIAIJ     *aij = (Mat_MPIAIJ*)mat->data;
        Mat            A = aij->A;
        Mat            B = aij->B;
        PetscInt       i,mA,nA, mB,nB;
        Vec            vA,vB;
        PetscScalar    *_vA, *_vB, *_v;
        PetscInt       *idxA,*idxB;
        PetscScalar    a,b;
        PetscReal      ta,tb;
	MatInfo        info;

        MatGetLocalSize( A, &mA, &nA );
        MatGetLocalSize( B, &mB, &nB );

        MatGetVecs( A, PETSC_NULL, &vA );
        MatGetVecs( B, PETSC_NULL, &vB );
        VecSet( vA, 0.0 );
        VecSet( vB, 0.0 );


        if (idx!=PETSC_NULL) {
                PetscMalloc( sizeof(PetscInt) * mA, &idxA );
                PetscMalloc( sizeof(PetscInt) * mB, &idxB );
        }
        else {
                idxA = idxB = PETSC_NULL;
        }

        MatGetInfo( A, MAT_GLOBAL_MAX, &info );
        if( info.nz_allocated > 0.0 ) {   MatGetRowMaxAbs( A, vA, idxA );   }

        MatGetInfo( B, MAT_GLOBAL_MAX, &info );
        if( info.nz_allocated > 0.0 ) {   MatGetRowMaxAbs( B, vB, idxB );   }


        VecGetArray( vA, &_vA );
        VecGetArray( vB, &_vB );
        VecGetArray( v, &_v );

        for (i=0; i<mA; i++) {
                a = _vA[i];
                b = _vB[i];

                ta = PetscAbsScalar(a);
                tb = PetscAbsScalar(b);

                /* compare  (ta < tb) */
                if (ta < tb) { //
                        _v[i] = _vA[i];
                        if (idx!=PETSC_NULL) {          idx[i] = idxA[i];               }
                } //
                else {
                        _v[i] = _vB[i];
                        if (idx!=PETSC_NULL) {          idx[i] = idxB[i];               }
                }
        }


        /* tidy up */
        if (idx!=PETSC_NULL) {
                PetscFree( idxA );
                PetscFree( idxB );
        }

        VecRestoreArray( vA, &_vA );
        VecRestoreArray( vB, &_vB );
        VecRestoreArray( v, &_v );

        Stg_VecDestroy( &vA );
        Stg_VecDestroy( &vB );

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMinAbs_MATAIJ(Mat A,Vec r,PetscInt idx[])
{
        PetscInt i,j;
        PetscInt m, n;

        const PetscInt *cols;
        PetscInt ncols,r_col;
        const PetscScalar *vals;
        PetscScalar *r_vals;
        PetscReal cmp, abv;


        MatGetOwnershipRange( A, &m, &n );

        VecGetArray( r, &r_vals );
        for( i=m; i<n; i++ ) {
                MatGetRow( A, i, &ncols, &cols, &vals );
                if( ncols==0 ) continue;

                cmp = PetscAbsScalar(vals[0]);
                r_col = cols[0];
                for( j=1; j<ncols; j++ ) {
                        abv = PetscAbsScalar(vals[j]);
                        if( abv < cmp ) {  /* Look for minimum */
                                cmp = abv;
                                r_col = cols[j];
                        }
                }
                if( idx ) {  idx[i-m] = r_col;  }
                r_vals[i-m] = cmp;

                MatRestoreRow( A, i, &ncols, &cols, &vals );
        }
        VecRestoreArray( r, &r_vals );

        PetscFunctionReturn(0);
}

PetscErrorCode PetscExtMatGetRowMinAbs(Mat mat,Vec vec,PetscInt idx[])
{
        PetscTruth is_mpiaij;

        PetscTypeCompare( (PetscObject)mat, MATMPIAIJ, &is_mpiaij );
        if (is_mpiaij==PETSC_TRUE) {    PetscExtMatGetRowMinAbs_MATAIJ(mat,vec,idx);       }
        else {                          PetscExtMatGetRowMinAbs_MATAIJ(mat,vec,idx);       }

        PetscFunctionReturn(0);
}

/* I DON' RECOMMEND USING THIS */
PetscErrorCode PetscExtLoadMPIAIJOperations(Mat mat)
{
	PetscTruth is_mpiaij;
	
	PetscTypeCompare( (PetscObject)mat, MATMPIAIJ, &is_mpiaij );
	if (is_mpiaij!=PETSC_TRUE) {
		PetscFunctionReturn(0);
	}
	
	mat->ops->getrowmax = PetscExtMatGetRowMax_MPIAIJ;
	mat->ops->getrowmin = PetscExtMatGetRowMin_MPIAIJ;
	mat->ops->getrowmaxabs = PetscExtMatGetRowMaxAbs_MPIAIJ;
	
	
	PetscFunctionReturn(0);

}
