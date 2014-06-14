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


/*
Implements a basic block vector implementation.

Allows the definition of a block vector of the form
{X}^T = { {x1}, {x2}, ..., {xn} }^T
where {x1} is a vector of any type.

It is possible some methods implemented will not support 
nested block vectors. That is {x1} cannot be of type "block".
More testing needs to be done to verify this.

*/

#include <stdlib.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <private/vecimpl.h>

#include "vec_block_impl.h"
#include "private/vec/petscvec-block.h"


#undef __FUNCT__  
#define __FUNCT__ "VecSetUp_Block"
PetscErrorCode VecSetUp_Block( Vec V )
{
	Vec_Block *ctx = (Vec_Block*)V->data;
	PetscInt i;
	MPI_Comm comm;
	
	if( ctx->setup_called == PETSC_TRUE ) PetscFunctionReturn(0);
	
	
	ctx->nb = V->map->N;
	if( ctx->nb < 0 ) {
		Stg_SETERRQ( PETSC_ERR_ARG_WRONG, "Cannot create VEC_BLOCK with < 0 blocks." );
	}
	
	/* Communicator for sub vectors is same as block vector */
	PetscObjectGetComm( (PetscObject)V, &comm );
	
	/* Create space */
	ctx->v = (Vec*)malloc( sizeof(Vec) * ctx->nb );
	
	for( i=0; i<ctx->nb; i++ ) {
		ctx->v[i] = PETSC_NULL;
		/* Do not allocate memory for internal sub blocks */
		//VecCreate( comm, &ctx->v[i] );
	}
	
	ctx->setup_called = PETSC_TRUE;
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "VecDestroy_Block"
PetscErrorCode VecDestroy_Block( Vec v )
{
	Vec_Block        *vs = (Vec_Block*)v->data;
	PetscErrorCode ierr;
	PetscInt i;
	
	PetscFunctionBegin;
	
	/* if memory was published with AMS then destroy it */
	ierr = PetscObjectDepublish(v);CHKERRQ(ierr);
	
	if( vs->v != PETSC_NULL ) {
		for( i=0; i<vs->nb; i++ ) {
			if( vs->v[i] != PETSC_NULL ) {
                Stg_VecDestroy( &(vs->v[i]) );
                vs->v[i] = PETSC_NULL;
			}
		}
		free( vs->v );
	}
	
	
	ierr=PetscFree( vs );CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecPublish_Block"
PetscErrorCode VecPublish_Block(PetscObject obj)
{
	PetscFunctionBegin;
	PetscFunctionReturn(0);
}




/* ================================== Private functions ===================================== */


#undef __FUNCT__  
#define __FUNCT__ "_VecBlockGetValues"
PetscErrorCode _VecBlockGetValues( 
		Vec x,
		PetscInt m, const PetscInt idxm[],
		Vec vec[] )
{
	Vec_Block *b = (Vec_Block*)x->data;
	PetscInt i;
	PetscInt row;
	
	if (!m ) return(0);
	
	for( i=0; i<m; i++ ) {
		row = idxm[i];
		if( row >= x->map->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,x->map->N-1);
		
		vec[ i ] = b->v[row]; // if( row_oriented )
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "_VecBlockGetValue"
PetscErrorCode _VecBlockGetValue( Vec x, PetscInt idxm, Vec *vec )
{
	_VecBlockGetValues( x, 1, &idxm, vec );
	PetscFunctionReturn(0);
}


/* ================================== Public functions ===================================== */

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockSetValues_Block"
PetscErrorCode VecBlockSetValues_Block( 
		Vec V, 
		PetscInt m, const PetscInt idxm[],
		const Vec vec[], InsertMode addv )
{
	Vec_Block *b = (Vec_Block*)V->data;;
	PetscInt i;
	PetscInt row;
	
	
	if (!m ) PetscFunctionReturn(0);
	
	for( i=0; i<m; i++ ) {
		row = idxm[i];
		if( row >= V->map->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,V->map->N-1);
		if( b->v[row] != PETSC_NULL ) {
           Stg_VecDestroy( &(b->v[row]) );
		}
	}
	
	for( i=0; i<m; i++ ) {
		row = idxm[i];
		if( row >= V->map->N ) Stg_SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,V->map->N-1);
		
		if( addv == INSERT_VALUES ) {
			b->v[row] = vec[ i ];
			PetscObjectReference( (PetscObject)vec[i] );
		}
		else {
			VecAXPY( b->v[row], 1.0, vec[ i] );
		}
	}
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecBlockSetValues"
PetscErrorCode VecBlockSetValues( 
		Vec V, 
		PetscInt m, const PetscInt idxm[],
		const Vec vec[], InsertMode addv )
{
	PetscErrorCode ierr,(*f)(Vec,PetscInt,const PetscInt*,const Vec*,InsertMode);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)V,"VecBlockSetValues_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(V,m,idxm,vec,addv);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}



/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockSetValue_Block"
PetscErrorCode VecBlockSetValue_Block( Vec x, PetscInt idxm, Vec vec, InsertMode addv )
{
	VecBlockSetValues( x, 1, &idxm, &vec, addv );
	PetscFunctionReturn(0);	
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecBlockSetValue"
PetscErrorCode VecBlockSetValue( Vec x, PetscInt idxm, Vec vec, InsertMode addv )
{
	PetscErrorCode ierr,(*f)(Vec,PetscInt,Vec,InsertMode);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)x,"VecBlockSetValue_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(x,idxm,vec,addv);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}



/*
The functions VecBlockGetValues() and VecBlockGetValue() should possibly be renamed
to VecBlockGetSubVectors() to be more consistent with VecGetArray() / VecRestoreArray()

If either
	i) VecBlockGetSubVector()
	ii) VecBlockGetSubVectors()
are called, you MUST call
	VecBlockRestoreSubVectors()
otherwise the block vectors state will not be increased, even through its 
subvectors (subordinates) have potentially been modified and have had their
state increased.
*/

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockGetSubVector_Block"
PetscErrorCode VecBlockGetSubVector_Block( Vec X, PetscInt idxm, Vec *sx )
{
	_VecBlockGetValue( X, idxm, sx );
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecBlockGetSubVector"
PetscErrorCode VecBlockGetSubVector( Vec X, PetscInt idxm, Vec *sx )
{
	PetscErrorCode ierr,(*f)(Vec,PetscInt,Vec*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)X,"VecBlockGetSubVector_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(X,idxm,sx);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockGetSubVectors_Block"
PetscErrorCode VecBlockGetSubVectors_Block( Vec X, Vec **sx )
{
	Vec_Block *b;
	
	b = (Vec_Block*)X->data;
	
	*sx = b->v;
	
	PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecBlockGetSubVectors"
PetscErrorCode PETSCMAT_DLLEXPORT VecBlockGetSubVectors( Vec x, Vec **sx )
{
	PetscErrorCode ierr,(*f)(Vec,Vec**);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)x,"VecBlockGetSubVectors_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(x,sx);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockRestoreSubVectors_Block"
PetscErrorCode VecBlockRestoreSubVectors_Block( Vec X )
{
	PetscObjectStateIncrease((PetscObject)X);
	PetscFunctionReturn(0);
}
EXTERN_C_END
#undef __FUNCT__  
#define __FUNCT__ "VecBlockRestoreSubVectors"
PetscErrorCode PETSCMAT_DLLEXPORT VecBlockRestoreSubVectors( Vec x )
{
	PetscErrorCode ierr,(*f)(Vec);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)x,"VecBlockRestoreSubVectors_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(x);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */

void _do_sum_entries( Vec x, PetscInt *n )
{
	PetscTruth isblock;
	Vec_Block *ctx;
	PetscInt i,N;
	
	
	if( x == PETSC_NULL ) return;
	
	PetscTypeCompare( (PetscObject)x, "block", &isblock );
	if( isblock == PETSC_TRUE ) {
		ctx = (Vec_Block*)x->data;
		if( ctx->v == PETSC_NULL ) return;
		
		for( i=0; i<ctx->nb; i++ ) {
			_do_sum_entries( ctx->v[i], n );
		}
	}
	else {
		VecGetSize( x, &N );
		*n = *n + N;
		return;
	}
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockSumEntries_Block"
PetscErrorCode VecBlockSumEntries_Block( Vec X, PetscInt *n )
{
	Vec_Block *ctx = (Vec_Block*)X->data;
	PetscInt i;
	
	*n = 0;
	for( i=0; i<ctx->nb; i++ ) {
		_do_sum_entries( ctx->v[i], n );
	}
	
	PetscFunctionReturn(0);
}
EXTERN_C_END
		
#undef __FUNCT__  
#define __FUNCT__ "VecBlockSumEntries"
PetscErrorCode VecBlockSumEntries( Vec X, PetscInt *n )
{
	PetscErrorCode ierr,(*f)(Vec,PetscInt*);
	
	PetscFunctionBegin;
	ierr = PetscObjectQueryFunction((PetscObject)X,"VecBlockSumEntries_C",(void (**)(void))&f);CHKERRQ(ierr);
	if (f) {
		ierr = (*f)(X,n);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecCreateFromVecBlock_Block"
PetscErrorCode VecCreateFromVecBlock_Block( Vec x, const VecType type, Vec *mX )
{
  Vec_Block *bx = (Vec_Block*)x->data;
  PetscTruth is_block;
  PetscInt i;
  PetscInt sum, len;
  MPI_Comm comm;      

  PetscTypeCompare( (PetscObject)x, "block", &is_block );
  if( !is_block ) Stg_SETERRQ( PETSC_ERR_SUP, "Cannot merge non-block vectors \n" ); 
        
       
  /* Get block size on all procs */
  sum = 0;
  for( i=0; i<bx->nb; i++ ) {
    /* check sub vectors are not blocks themselves */
    PetscTypeCompare( (PetscObject)bx->v[i], "block", &is_block );
    if( is_block ) Stg_SETERRQ( PETSC_ERR_SUP, "Cannot merge sub vectors which are themselves block vectors \n" );
    VecGetSize( bx->v[i], &len );
    sum = sum + len;
  }
        
  PetscObjectGetComm( (PetscObject)x, &comm );
  VecCreate( comm, mX );
  VecSetSizes( *mX, PETSC_DECIDE, sum );
  if( type!=PETSC_NULL ) {
    VecSetType( *mX, type );
  }
  VecSetFromOptions(*mX);    
    
  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecCreateFromVecBlock"
PetscErrorCode VecCreateFromVecBlock( Vec x, const VecType type, Vec *mX )
{
        PetscErrorCode ierr,(*f)(Vec,const VecType,Vec*);

        PetscFunctionBegin;
        ierr = PetscObjectQueryFunction((PetscObject)x,"VecCreateFromVecBlock_C",(void (**)(void))&f);CHKERRQ(ierr);
        if (f) {
                ierr = (*f)(x,type,mX);CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
}

/* =========================================================================== */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecBlockMergeSubVectors_Block"
PetscErrorCode VecBlockMergeSubVectors_Block( Vec x, InsertMode addv, const VecType type, Vec *mX )
{
  Vec_Block *bx = (Vec_Block*)x->data;
  VecScatter *vscat;
  int *block_offset;
  PetscInt i,local_len,sum,Ni;
  IS is; 
  Vec merged;
  MPI_Comm comm;


  /* check if mX exists */
  if( *mX == PETSC_NULL ) {
    VecCreateFromVecBlock( x, type, mX );
  }
  else {
    PetscInt N_x, N_mX;
    PetscTruth is_same;

    /* check type matches existing vector */
    is_same = PETSC_FALSE;
    PetscTypeCompare( (PetscObject)(*mX), type, &is_same );
    if(is_same==PETSC_FALSE) {
      Stg_SETERRQ( PETSC_ERR_ARG_NOTSAMETYPE, "Vector being re-used is of wrong type");
    }

    /* check length is consistent */
    VecBlockSumEntries( x, &N_x );
    VecGetSize( *mX, &N_mX );
    if( N_x!=N_mX ) { Stg_SETERRQ( PETSC_ERR_ARG_SIZ, "Vector being re-used is wrong global length"); }

  }
  merged = *mX; 

  /* get global offests for the merged vector */        
  PetscMalloc( sizeof(int)*bx->nb, &block_offset );
  sum = 0;
  for( i=0; i<bx->nb; i++ ) {
    VecGetSize( bx->v[i], &Ni );
    block_offset[i] = sum;
    sum = sum + Ni;
  }

  PetscObjectGetComm( (PetscObject)x, &comm );

  /* generate the scatter */
  PetscMalloc( sizeof(VecScatter)*bx->nb, &vscat );                 
             
  for( i=0; i<bx->nb; i++ ) {
    PetscInt s,e;

    VecGetLocalSize( bx->v[i], &local_len );
    VecGetOwnershipRange( bx->v[i], &s, &e );
                        
    ISCreateStride( comm, local_len, block_offset[i]+s, 1, &is );
    VecScatterCreate( bx->v[i], PETSC_NULL, merged, is, &vscat[i] ); 
    Stg_ISDestroy( &is );
  }

  /* do the scatter */
  for( i=0; i<bx->nb; i++ ) {
    VecScatterBegin( vscat[i], bx->v[i], merged, addv, SCATTER_FORWARD );
    VecScatterEnd( vscat[i], bx->v[i], merged, addv, SCATTER_FORWARD );
  }

  /* tidy up */
  PetscFree( block_offset );

  for( i=0; i<bx->nb; i++ ) {
    Stg_VecScatterDestroy( &(vscat[i]) );
  }
  PetscFree( vscat );


  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "VecBlockMergeSubVectors"
PetscErrorCode VecBlockMergeSubVectors( Vec x, InsertMode addv, const VecType type, Vec *mX )
{
        PetscErrorCode ierr,(*f)(Vec,InsertMode,const VecType,Vec*);

        PetscFunctionBegin;
        ierr = PetscObjectQueryFunction((PetscObject)x,"VecBlockMergeSubVectors_C",(void (**)(void))&f);CHKERRQ(ierr);
        if (f) {
                ierr = (*f)(x,addv,type,mX);CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
}


/* Constructors */
PetscErrorCode VecBlockSetOperations( struct _VecOps *o )
{

/* 0 */     o->duplicate = VecDuplicate_Block;
                        o->duplicatevecs=VecDuplicateVecs_Default;
                        o->destroyvecs=VecDestroyVecs_Default;
                        o->dot=VecDot_Block;
                        o->mdot = VecMDot_Block;

/* 5 */         o->norm = VecNorm_Block;
                        o->tdot=VecTDot_Block;
                        o->mtdot = VecMTDot_Block;
                        o->scale=VecScale_Block;
                        o->copy=VecCopy_Block;

/* 10 */        o->set=VecSet_Block;
                        o->swap=VecSwap_Block;
                        o->axpy=VecAXPY_Block;
                        o->axpby=VecAXPBY_Block;
                        o->maxpy=VecMAXPY_Block;

/* 15 */        o->aypx=VecAYPX_Block;
                        o->waxpy=VecWAXPY_Block;
  o->axpbypcz = PETSC_NULL;
                        o->pointwisemult = VecPointwiseMult_Block;
                        o->pointwisedivide=VecPointwiseDivide_Block;
/* 20 */ o->setvalues=__VecSetValues_Block;
        o->assemblybegin = 0; //VecAssemblyBegin_Empty;
                        o->assemblyend=0; //VecAssemblyEnd_Empty;
                        o->getarray=__VecGetArray_Block;
                        o->getsize=VecGetSize_Block; //VecGetSize_Empty;

/* 25 */                        o->getlocalsize=VecGetSize_Block; //VecGetLocalSize_Empty;
       o->restorearray= __VecRestoreArray_Block;
                        o->max=VecMax_Block;
                        o->min=VecMin_Block;
                        o->setrandom=__VecSetRandom_Block;

/* 30 */                        o->setoption=0; //VecSetOption_Empty;
        o->setvaluesblocked=__VecSetValuesBlocked_Block;
                        o->destroy=VecDestroy_Block;
                        o->view=VecView_Block;
                        o->placearray = __VecPlaceArray_Block;

/* 35 */                        o->replacearray=__VecReplaceArray_Block;
        o->dot_local=VecDot_Block; //VecDotLocal_Empty;
                        o->tdot_local=VecTDot_Block; //VecTDotLocal_Empty;
                        o->norm_local=VecNorm_Block; //VecNormLocal_Empty;
                        o->mdot_local=VecMDot_Block; //VecMDotLocal_Empty;

/* 40 */                        o->mtdot_local=VecMTDot_Block; //VecMTDotLocal_Empty;
        o->loadintovector=__VecLoadIntoVector_Block;
                        o->reciprocal=VecReciprocal_Block;
                        //o->viewnative=0; //VecViewNative_Empty;
                        o->conjugate=VecConjugate_Block;

/* 45 */                        o->setlocaltoglobalmapping=0; //VecSetLocalToGlobalMapping_Empty;
  o->setvalueslocal=0; //VecSetValuesLocal_Empty;
                       o->resetarray= __VecResetArray_Block;
                        o->setfromoptions=0; //VecSetFromOptions_Empty;
                        o->maxpointwisedivide=VecMaxPointwiseDivide_Block;

/* 50 */                        o->load=0; //VecLoad_Empty;
        o->pointwisemax=__VecPointwiseMax_Block;
                        o->pointwisemaxabs=__VecPointwiseMaxAbs_Block;
                        o->pointwisemin=__VecPointwiseMin_Block;
                        o->getvalues=__VecGetValues_Block;



  PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecCreate_Block"
PetscErrorCode VecCreate_Block( Vec V )
{
	Vec_Block      *s;
	PetscErrorCode ierr;
	
	
	PetscFunctionBegin;
	
	
//	ierr = PetscMemcpy( V->ops, &VecBlockOps_Values,sizeof(struct _VecOps) );CHKERRQ(ierr);
	ierr = VecBlockSetOperations( V->ops );

//	V->precision       = PETSC_SCALAR;
//	V->bops->publish   = VecPublish_Block;
	V->petscnative     = PETSC_TRUE;
	
	/* allocate and set pointer for implememtation data */
	ierr = PetscMalloc( sizeof(Vec_Block), &s );CHKERRQ(ierr);
	ierr = PetscMemzero( s, sizeof(Vec_Block) );CHKERRQ(ierr);
	V->data            = (void*)s;

	ierr = PetscObjectChangeTypeName((PetscObject)V,"block");CHKERRQ(ierr);
	
	s->setup_called  = PETSC_FALSE;
	s->nb            = -1;
	s->v             = PETSC_NULL;
	
	/*
	M = V->map.N;
	if (V->map.bs == -1) V->map.bs = 1;
	PetscMapSetUp( &V->map );
	PetscMapSetLocalSize( &V->map, M );
	PetscMapSetSize( &V->map, M );
	*/
	
	VecSetUp_Block( V );
	
	ierr = PetscPublishAll(V);CHKERRQ(ierr);
	
	
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockRestoreSubVectors_C",
			"VecBlockRestoreSubVectors_Block",
			VecBlockRestoreSubVectors_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockGetSubVectors_C",
			"VecBlockGetSubVectors_Block",
			VecBlockGetSubVectors_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockGetSubVector_C",
			"VecBlockGetSubVector_Block",
			VecBlockGetSubVector_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockSetValue_C",
			"VecBlockSetValue_Block",
			VecBlockSetValue_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockSetValues_C",
			"VecBlockSetValues_Block",
			VecBlockSetValues_Block);CHKERRQ(ierr);
	
	ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockSumEntries_C",
			"VecBlockSumEntries_Block",
			VecBlockSumEntries_Block);CHKERRQ(ierr);
	
        ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecCreateFromVecBlock_C",
                        "VecCreateFromVecBlock_Block",
                        VecCreateFromVecBlock_Block);CHKERRQ(ierr);
        ierr = PetscObjectComposeFunctionDynamic((PetscObject)V,"VecBlockMergeSubVectors_C",
                        "VecBlockMergeSubVectors_Block",
                        VecBlockMergeSubVectors_Block );CHKERRQ(ierr);

	
	PetscFunctionReturn(0);
}
EXTERN_C_END

