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

#include <petsc.h>
#include <petscvec.h>
#include <petscversion.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >=3) )
  #include <petsc-private/vecimpl.h>
#else
  #include <private/vecimpl.h>
#endif

#include "vec_block_impl.h"
#include "private/vec/petscvec-block.h"




PetscErrorCode BlockCheck_2( Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	
	if( bx->setup_called == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vector x not setup." );

	}
	if( by->setup_called == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vector y not setup." );
	}
	if( bx->nb != by->nb ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vectors have different numbers of blocks." );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode BlockCheck_3( Vec w, Vec x, Vec y )
{
	Vec_Block *bw = (Vec_Block*)w->data;
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	
	if( bw->setup_called == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vector w not setup." );
	}
	if( bx->setup_called == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vector x not setup." );
	}
	if( by->setup_called == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vector y not setup." );
	}
	
	
	
	if( bx->nb != by->nb ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vectors x and y have different numbers of blocks." );
	}
	if( bx->nb != bw->nb ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vectors x and w have different numbers of blocks." );
	}
	if( bw->nb != by->nb ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Block vectors w and y have different numbers of blocks." );
	}
	
	PetscFunctionReturn(0);
}


/* supports nested blocks */
PetscErrorCode VecCopy_Block( Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i;
	
	
	BlockCheck_2(x,y);
	for( i=0; i<bx->nb; i++ ) {
		VecCopy( bx->v[i], by->v[i] );
	}
	
	PetscFunctionReturn(0);
}


/* supports nested blocks */
PetscErrorCode VecDuplicate_Block( Vec x, Vec *y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	//Vec_Block *by;
	PetscInt i;
	Vec _y;
	Vec Y;
	
	
	VecCreate( ((PetscObject)x)->comm, &_y );
	VecSetSizes( _y, bx->nb, bx->nb );
	VecSetType( _y, "block" );
	VecSetUp( _y );

	//by = (Vec_Block*)_y->data;
	for( i=0; i<bx->nb; i++ ) {
		VecDuplicate( bx->v[i], &Y );
		/* Only required to Stg_VecDestroy( &) IF block allocates memory, should use the BlockSetValue to be safe */
	//	Stg_VecDestroy( & by->v[i] );
	//	by->v[i] = Y;
		VecBlockSetValue( _y, i, Y, INSERT_VALUES );
		Stg_VecDestroy( &Y ); /* Hand over control of Y to the block vector _y */
	}
	
	*y = _y;
	
	PetscFunctionReturn(0);
}



/* supports nested blocks */
PetscErrorCode VecDot_Block( Vec x, Vec y, PetscScalar *val )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	PetscScalar	x_dot_y, _val;
	
	nr = bx->nb;
	_val = 0.0;
	for( i=0; i<nr; i++ ) {
		VecDot( bx->v[i], by->v[i], &x_dot_y );
		_val = _val + x_dot_y;
	}
	(*val) = _val;
	PetscFunctionReturn(0);
}

/* supports nested blocks */
PetscErrorCode VecTDot_Block( Vec x, Vec y, PetscScalar *val )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	PetscScalar	x_dot_y, _val;
	
	nr = bx->nb;
	_val = 0.0;
	for( i=0; i<nr; i++ ) {
		VecTDot( bx->v[i], by->v[i], &x_dot_y );
		_val = _val + x_dot_y;
	}
	(*val) = _val;
	PetscFunctionReturn(0);
}



PetscErrorCode VecAXPY_Block( Vec y, PetscScalar alpha, Vec x )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecAXPY( by->v[i], alpha, bx->v[i] );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecAYPX_Block( Vec y, PetscScalar alpha, Vec x )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecAYPX( by->v[i], alpha, bx->v[i] );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecAXPBY_Block( Vec y, PetscScalar alpha, PetscScalar beta, Vec x )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecAXPBY( by->v[i], alpha, beta, bx->v[i] );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecScale_Block( Vec x, PetscScalar alpha )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecScale( bx->v[i], alpha );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecPointwiseMult_Block( Vec w, Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	Vec_Block *bw = (Vec_Block*)w->data;
	PetscInt i, nr;
	
	
	BlockCheck_3(w,x,y);
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecPointwiseMult( bw->v[i], bx->v[i], by->v[i] );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecPointwiseDivide_Block( Vec w, Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	Vec_Block *bw = (Vec_Block*)w->data;
	PetscInt i, nr;
	
	
	BlockCheck_3(w,x,y);
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecPointwiseDivide( bw->v[i], bx->v[i], by->v[i] );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecReciprocal_Block( Vec x )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecReciprocal( bx->v[i] );
	}
	
	PetscFunctionReturn(0);
}


PetscErrorCode VecNorm_Block( Vec xin, NormType type, PetscReal* z )
{
	Vec_Block *bx = (Vec_Block*)xin->data;
	PetscInt i, nr;
	PetscReal z_i;
	PetscReal _z;
	
	nr = bx->nb;
	_z = 0.0;
	
/*	
	for( i=0; i<nr; i++ ) {
		VecNorm( bx->v[i], type, &z_i );
		
		if( type == NORM_2 ) {
			_z = _z + z_i*z_i;
		}
		else if( type == NORM_1 ) {
			_z = _z + z_i;
		}
		else if( type == NORM_INFINITY ) {
			if( z_i > _z ) {
				_z = z_i;
			}
		}
	}
	if( type == NORM_2 ) {
		_z = sqrt(_z);
	}
*/
	
	if( type == NORM_2 ) {
		PetscScalar dot;
#ifdef PETSC_USE_COMPLEX
		PetscReal im,re;
#endif
		
		VecDot( xin, xin, &dot );
#ifdef PETSC_USE_COMPLEX
		re = PetscRealPart( dot );
		im = PetscImaginaryPart( dot );
		_z = sqrt( re - im );
#else
		_z = sqrt(dot);
#endif 
	}
	else if( type == NORM_1 ) {
		for( i=0; i<nr; i++ ) {
			VecNorm( bx->v[i], type, &z_i );
			_z = _z + z_i;
		}
	}
	else if( type == NORM_INFINITY ) {
		for( i=0; i<nr; i++ ) {
			VecNorm( bx->v[i], type, &z_i );
			if( z_i > _z ) {
				_z = z_i;
			}
		}
	}
	
	*z = _z;
	
	PetscFunctionReturn(0);
}


PetscErrorCode VecMAXPY_Block( Vec y, PetscInt nv, const PetscScalar alpha[], Vec *x )
{
	PetscInt v;
	
	for( v=0; v<nv; v++ ) {
		/* Do axpy on each vector, v */
		VecAXPY( y, alpha[v], x[v] );
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecMDot_Block( Vec x, PetscInt nv, const Vec y[], PetscScalar *val )
{
	PetscInt j;
	
	for( j=0; j<nv; j++ ) {
		VecDot( x, y[j], &val[j] );
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecMTDot_Block( Vec x, PetscInt nv, const Vec y[], PetscScalar *val )
{
	PetscInt j;
	
	for( j=0; j<nv; j++ ) {
		VecTDot( x, y[j], &val[j] );
	}
	
	PetscFunctionReturn(0);
}


PetscErrorCode VecSet_Block( Vec x, PetscScalar alpha )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	PetscInt i, nr;
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecSet( bx->v[i], alpha );
	}
	PetscFunctionReturn(0);
}

PetscErrorCode VecConjugate_Block( Vec x )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	PetscInt j, nr;
	
	nr = bx->nb;
	for( j=0; j<nr; j++ ) {
		VecConjugate( bx->v[j] );
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecSwap_Block( Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	
	
	BlockCheck_2(x,y);
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecSwap( bx->v[i], by->v[i] );
	}
	
	PetscFunctionReturn(0);
}

PetscErrorCode VecWAXPY_Block( Vec w, PetscScalar alpha, Vec x, Vec y )
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	Vec_Block *bw = (Vec_Block*)w->data;
	PetscInt i, nr;
	
	
	BlockCheck_3(w,x,y);
	
	nr = bx->nb;
	for( i=0; i<nr; i++ ) {
		VecWAXPY( bw->v[i], alpha, bx->v[i], by->v[i] );
	}
	PetscFunctionReturn(0);
}

void __vec_max_block( Vec x, PetscInt *cnt, PetscInt *p, PetscReal *max )
{
	Vec_Block   *bx = (Vec_Block*)x->data;

	PetscInt     i, nr;
	PetscTruth   isblock;
	PetscInt     L;
	PetscInt     _entry_loc;
	PetscReal    _entry_val;
	
	
	Stg_PetscTypeCompare( (PetscObject)x, "block", &isblock );
	if( isblock == PETSC_FALSE ) {
		/* Not block */
		VecMax( x, &_entry_loc, &_entry_val );
		if( _entry_val > *max ) {
			*max = _entry_val;
			*p = _entry_loc + *cnt;
		}
		
		VecGetSize( x, &L );
		*cnt = *cnt + L;
		
		return;
	}
	
	/* Otherwise we have a block */
	bx = (Vec_Block*)x->data;
	nr = bx->nb;
	
	/* now descend recursively */
	for( i=0; i<nr; i++ ) {
		__vec_max_block( bx->v[i], cnt, p, max );
	}
	
	return;
}


/* supports nested blocks */
PetscErrorCode VecMax_Block( Vec x, PetscInt *p, PetscReal *max )
{
	PetscInt  cnt;
	

	cnt = 0;
	*p = 0;
	*max = 0.0;
	__vec_max_block( x, &cnt, p, max );
	
	PetscFunctionReturn(0);
}



void __vec_min_block( Vec x, PetscInt *cnt, PetscInt *p, PetscReal *min )
{
	Vec_Block   *bx = (Vec_Block*)x->data;

	PetscInt     i, nr;
	PetscTruth   isblock;
	PetscInt     L;
	PetscInt     _entry_loc;
	PetscReal    _entry_val;
	
	
	Stg_PetscTypeCompare( (PetscObject)x, "block", &isblock );
	if( isblock == PETSC_FALSE ) {
		/* Not block */
		VecMin( x, &_entry_loc, &_entry_val );
		if( _entry_val < *min ) {
			*min = _entry_val;
			*p = _entry_loc + *cnt;
		}
		
		VecGetSize( x, &L );
		*cnt = *cnt + L;
		
		return;
	}
	
	/* Otherwise we have a block */
	bx = (Vec_Block*)x->data;
	nr = bx->nb;
	
	/* now descend recursively */
	for( i=0; i<nr; i++ ) {
		__vec_min_block( bx->v[i], cnt, p, min );
	}
	
	return;
}

PetscErrorCode VecMin_Block( Vec x, PetscInt *p, PetscReal *min )
{
	PetscInt  cnt;
	
	cnt = 0;
	*p = 0;
	*min = 1.0e308;
	__vec_min_block( x, &cnt, p, min );
	
	PetscFunctionReturn(0);
}


/* supports nested blocks */
void print_vec_structure( Vec x, int parent_index, PetscViewer viewer )
{
	Vec_Block	*bx;
	VecType		type;
	PetscInt	m;
	int i;
	const char *name;
	PetscTruth is_block;
	
	
	VecGetType( x, &type );
	name = ((PetscObject)x)->prefix;
	
	
	is_block = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)x, "block", &is_block );
	if( is_block  == PETSC_FALSE ) {
		VecGetSize( x, &m );
		
		PetscViewerASCIIPushTab( viewer );
		
		if( name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%p), type=%s, rows=%d \n", parent_index,0, x, type, m ); }
		else { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%s), type=%s, rows=%d \n", parent_index,0, name, type, m ); }
		
		PetscViewerASCIIPopTab( viewer );
		return;
	}
	else {
		VecGetSize( x, &m );
		
		PetscViewerASCIIPushTab( viewer );
		
//		PetscViewerASCIIPrintf( viewer, "BVec(%p)- [%d] x [%d] \n", x,m, 1 );
		if( name == PETSC_NULL ) { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%p), type=%s, rows=%d \n", parent_index,0, x, type, m ); }
		else { PetscViewerASCIIPrintf( viewer, "(%d,%d) - (%s), type=%s, rows=%d \n", parent_index,0, name, type, m ); }
		
		bx = (Vec_Block*)x->data;
		for( i=0; i<bx->nb; i++ ) {
			print_vec_structure( bx->v[i],i, viewer );
		}
		PetscViewerASCIIPopTab( viewer );
	}
}

void print_vec_contents( Vec x, Vec parent_vec, PetscInt pid, PetscViewer viewer )
{
	Vec_Block	*bx;
	VecType		type;
	PetscInt	m;
	int i;
	const char *name, *pname;
	PetscTruth is_block;
	
	
	name = ((PetscObject)x)->prefix;
	if( parent_vec != PETSC_NULL ) {
		pname = ((PetscObject)parent_vec)->prefix;
	}
	else pname = PETSC_NULL;
	
	VecGetType( x, &type );
	is_block = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)x, "block", &is_block );
	
	if( is_block  == PETSC_FALSE ) {
		VecGetSize( x, &m );
		
		
		PetscViewerASCIIPushTab( viewer );
		
		//PetscViewerASCIIPrintf( viewer, "(i,1) - [sub vector]: (Name / ptr) , [parent vector]: (Name / ptr) \n" );
		PetscViewerASCIIPrintf( viewer, "-------------------- Begin ( %p ) --------------------\n", x );
		if( name == PETSC_NULL ) { 
			if( pname == PETSC_NULL ) PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub vec]: (%p), [parent vec]: %p \n", pid,0, x, parent_vec ); 
			else PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub vec]: (%p), [parent vec]: %s \n", pid,0, x, pname ); 
		}
		else { 
			if( pname == PETSC_NULL ) PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub vec]: (%s), [parent vec]: %p \n", pid,0, name, parent_vec ); 
			else PetscViewerASCIIPrintf( viewer, "(%d,%d) - [sub vec]: (%s), [parent vec]: %s \n", pid,0, name, pname ); 
		}
		VecView( x, viewer);
		PetscViewerASCIIPrintf( viewer, "--------------------   End ( %p ) --------------------\n\n", x );
		
		PetscViewerASCIIPopTab( viewer );
		
		return;
	}
	else {
		VecGetSize( x, &m );
		
		bx = (Vec_Block*)x->data;
		for( i=0; i<bx->nb; i++ ) {
			if( bx->v[i] == PETSC_NULL ) {
			}
			else {
				print_vec_contents( bx->v[i],x, i, viewer );
			}
		}
		
	}
}

/* supports nested blocks */
PetscErrorCode VecView_Block( Vec x, PetscViewer viewer )
{
	PetscViewerFormat format;
	PetscErrorCode ierr;
	Vec_Block *bx = (Vec_Block*)x->data;
	PetscTruth isascii;
		
	Stg_PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {

		PetscViewerASCIIPrintf( viewer, "Vector Object:\n");
		PetscViewerASCIIPushTab( viewer );		// push0
		PetscViewerASCIIPrintf( viewer, "type=block, rows=%d \n",bx->nb);
		
		PetscViewerASCIIPrintf(viewer,"VecBlock structure: \n" );
		PetscViewerASCIIPrintf( viewer, "(i,0) - (Name / ptr), type=, rows=\n" );
		print_vec_structure( x, 0, viewer );
		
		
		ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
		if( format == PETSC_VIEWER_ASCII_INFO_DETAIL ) {
			PetscViewerASCIIPrintf(viewer,"VecBlock contents: \n" );
			PetscViewerASCIIPrintf( viewer, "(i,0) - [sub vec]: (Name / ptr) , [parent vec]: (Name / ptr) \n" );
			PetscViewerASCIIPushTab( viewer );
			print_vec_contents( x, PETSC_NULL, 0, viewer );
			PetscViewerASCIIPopTab( viewer );
		}
	
	
		PetscViewerASCIIPopTab( viewer );		// pop0
	}
	PetscFunctionReturn(0);
}

/* Returns the number of blocks in size */
PetscErrorCode  VecGetSize_Block(Vec x,PetscInt *size)
{
	Vec_Block *bx = (Vec_Block*)x->data;
	
	*size = bx->nb;
	
	PetscFunctionReturn(0);
}



/* NOT SUPPORTED */
PetscErrorCode __VecSetValues_Block( Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora )
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}


PetscErrorCode __VecGetArray_Block(Vec x,PetscScalar *a[])
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}


PetscErrorCode __VecRestoreArray_Block(Vec x,PetscScalar *a[])
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}


PetscErrorCode  __VecSetRandom_Block(Vec x,PetscRandom rctx) 
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode  __VecSetValuesBlocked_Block(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode __VecPlaceArray_Block(Vec vec,const PetscScalar array[])
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode __VecReplaceArray_Block(Vec vec,const PetscScalar array[])
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode __VecLoadIntoVector_Block(PetscViewer viewer,Vec vec)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode __VecResetArray_Block(Vec vec)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode VecMaxPointwiseDivide_Block(Vec x,Vec y,PetscReal *max)
{
	Vec_Block *bx = (Vec_Block*)x->data;
	Vec_Block *by = (Vec_Block*)y->data;
	PetscInt i, nr;
	PetscReal local_max, m;
	
	BlockCheck_2( x,y );
	
	nr = bx->nb;
	m = 0.0;
	for( i=0; i<nr; i++ ) {
		VecMaxPointwiseDivide( bx->v[i], by->v[i], &local_max );
		if( local_max > m ) {
			m = local_max;
		}
	}
	(*max) = m;
	
	PetscFunctionReturn(0);
}


PetscErrorCode  __VecPointwiseMax_Block(Vec w,Vec x,Vec y)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

PetscErrorCode  __VecPointwiseMaxAbs_Block(Vec w,Vec x,Vec y)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}


PetscErrorCode  __VecPointwiseMin_Block(Vec w,Vec x,Vec y)
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}


PetscErrorCode  __VecGetValues_Block(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])
{
	Stg_SETERRQ( PETSC_ERR_SUP, "Operation not supported by VECBLOCK" );
	PetscFunctionReturn(0);
}

