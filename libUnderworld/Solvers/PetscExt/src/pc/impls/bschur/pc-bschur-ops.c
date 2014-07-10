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

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/kspimpl.h>
#else
  #include <private/kspimpl.h>
#endif
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/pcimpl.h>
#else
  #include <private/pcimpl.h>
#endif
#if ( (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3) )
  #include <petsc-private/matimpl.h>
#else
  #include <private/matimpl.h>
#endif

#include "private/vec/petscvec-block.h"
#include "private/mat/petscmat-block.h"
#include "private/mat/mat-schur.h"

#include "private/pc/pc-bschur.h"
#include "pc-bschur-impl.h"

/* externs */
extern const char *PCBSchurTypes[];




/* 0 */
PetscErrorCode PCSetUp_BSchur(PC pc)
{
	PC_BSchur   s = (PC_BSchur)pc->data;
	Mat         A, B;
	PetscTruth  is_block;
	PetscInt    M,N;
	Mat A11, A12, A21, A22;
	KSP ksp;
	
	
	
	PetscFunctionBegin; 
	
	
	/* check its okay to do stuff with this pc */
	Stg_PCGetOperators( pc, &A, &B, 0 );
	
	is_block = PETSC_FALSE;
	Stg_PetscTypeCompare( (PetscObject)B, "block", &is_block );
	if( is_block == PETSC_FALSE ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Preconditioner matrix B must be of type MAT_BLOCK" );
	}
	
	MatGetSize( B, &M, &N );
	if( (M != 2) || (N!=2) ) {
		Stg_SETERRQ( PETSC_ERR_SUP, "Preconditioner matrix B must be a 2 x 2 block matrix" );
	}
	
	MatBlockGetSubMatrix( B, 0,0, &A11 );
	MatBlockGetSubMatrix( B, 0,1, &A12 );
	MatBlockGetSubMatrix( B, 1,0, &A21 );
	MatBlockGetSubMatrix( B, 1,1, &A22 );
	
	
	
	
	if( s->explicit_operator == PETSC_TRUE ) {
		/* check we have enough matrices for application */
		if( s->application_type == PC_BSCHUR_UPPER ) {
			if( A11 == PETSC_NULL ) {
				Stg_SETERRQ( PETSC_ERR_USER, "You must provide A11 to define explicit matrix S_2 = A21 [diag(A11)]^{-1} A12");
			}
		}
		else if ( s->application_type == PC_BSCHUR_LOWER ) {
			if( A22 == PETSC_NULL ) {
				Stg_SETERRQ( PETSC_ERR_USER, "You must provide A22 to define explicit matrix S_1 = A12 [diag(A22)]^{-1} A21");
			}
		}
		
		
		
		
		if( s->explicit_schur == PETSC_NULL ) {
			
			if( s->application_type == PC_BSCHUR_UPPER ) {
				Vec diag;
				
				MatGetVecs( A11, &diag, 0 );
				MatGetDiagonal( A11, diag );
				VecReciprocal( diag );
				
			//	MatDuplicate( A12, &sA12 );
				MatDiagonalScale( A12, diag, 0 );
				
				MatMatMult( A21, A12, MAT_INITIAL_MATRIX, 1.0, &s->explicit_schur );
				MatAXPY( s->explicit_schur, -1.0, A22, DIFFERENT_NONZERO_PATTERN ); /* s = -A22 + s */
				
			//	Stg_MatDestroy( sA12 );
				VecReciprocal( diag );
				MatDiagonalScale( A12, diag, 0 );
				
				
				Stg_VecDestroy( &diag);
			}
			else if ( s->application_type == PC_BSCHUR_LOWER ) {
				Vec diag;
				
				MatGetVecs( A22, &diag, 0 );
				MatGetDiagonal( A22, diag );
				VecReciprocal( diag );
				
				MatDiagonalScale( A21, diag, 0 );
				
				MatMatMult( A12, A21, MAT_INITIAL_MATRIX, 1.0, &s->explicit_schur );
				MatAXPY( s->explicit_schur, -1.0, A11, DIFFERENT_NONZERO_PATTERN ); /* s = -A11 + s */
				
				VecReciprocal( diag );
				MatDiagonalScale( A21, diag, 0 );
				
				
				Stg_VecDestroy( &diag);
			}
			
		}
		else {
			/* === */
			if( s->application_type == PC_BSCHUR_UPPER ) {
				Vec diag;
				
				MatGetVecs( A11, &diag, 0 );
				MatGetDiagonal( A11, diag );
				VecReciprocal( diag );
				
			//	MatDuplicate( A12, &sA12 );
				MatDiagonalScale( A12, diag, 0 );
				
				MatMatMult( A21, A12, MAT_REUSE_MATRIX, 1.0, &s->explicit_schur );
				MatAXPY( s->explicit_schur, -1.0, A22, DIFFERENT_NONZERO_PATTERN ); /* s = -A22 + s */
				
			//	Stg_MatDestroy( sA12 );
				VecReciprocal( diag );
				MatDiagonalScale( A12, diag, 0 );
				
				Stg_VecDestroy( &diag);
			}
			else if ( s->application_type == PC_BSCHUR_LOWER ) {
				Vec diag;
				
				MatGetVecs( A22, &diag, 0 );
				MatGetDiagonal( A22, diag );
				VecReciprocal( diag );
				
				MatDiagonalScale( A21, diag, 0 );
				
				MatMatMult( A12, A21, MAT_INITIAL_MATRIX, 1.0, &s->explicit_schur );
				MatAXPY( s->explicit_schur, -1.0, A11, DIFFERENT_NONZERO_PATTERN ); /* s = -A11 + s */
				
				VecReciprocal( diag );
				MatDiagonalScale( A21, diag, 0 );
				
				
				Stg_VecDestroy( &diag);
			}
		}
		
		/* update operators for the block preconditioner */
		if( s->application_type == PC_BSCHUR_UPPER ) {
			Stg_KSPSetOperators( s->ksp_1, A11, A11, SAME_NONZERO_PATTERN );
			Stg_KSPSetOperators( s->ksp_2, s->explicit_schur, s->explicit_schur, SAME_NONZERO_PATTERN );
		}
		else if ( s->application_type == PC_BSCHUR_LOWER ) {
			Stg_KSPSetOperators( s->ksp_1, s->explicit_schur, s->explicit_schur, SAME_NONZERO_PATTERN );
			Stg_KSPSetOperators( s->ksp_2, A22, A22, SAME_NONZERO_PATTERN );
		}
	}
	else {
		
		if( s->schur == PETSC_NULL ) {
			const char *pc_prefix;
			
			if( s->application_type == PC_BSCHUR_UPPER ) {
				MatCreateSchur( PETSC_COMM_WORLD, A11,A12,A21,A22, 0, "MatSchur_A11", &s->schur );
			}
			else if ( s->application_type == PC_BSCHUR_LOWER ) {
				MatCreateSchur( PETSC_COMM_WORLD, A11,A12,A21,A22, 0, "MatSchur_A22", &s->schur );
			}
			
			MatSchurGetKSP( s->schur, &ksp );
			PCGetOptionsPrefix( pc, &pc_prefix );
			KSPSetOptionsPrefix( ksp, pc_prefix );
			KSPAppendOptionsPrefix( ksp, "pc_bschur_mat_schur_" );
		}
		
		/* operator operators for A11/A22 */
		MatSchurGetKSP( s->schur, &ksp );
		if( s->application_type == PC_BSCHUR_UPPER ) {
			Stg_KSPSetOperators( ksp, A11, A11, SAME_NONZERO_PATTERN );
		}
		else if ( s->application_type == PC_BSCHUR_LOWER ) {
			Stg_KSPSetOperators( ksp, A22, A22, SAME_NONZERO_PATTERN );
		}
		MatAssemblyBegin(s->schur, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd(s->schur, MAT_FINAL_ASSEMBLY );
		
		
		/* update the operators for the preconditioner */
		if( s->application_type == PC_BSCHUR_UPPER ) {
			Stg_KSPSetOperators( s->ksp_1, A11, A11, SAME_NONZERO_PATTERN );
			Stg_KSPSetOperators( s->ksp_2, s->schur, s->schur, SAME_NONZERO_PATTERN );
		}
		else if ( s->application_type == PC_BSCHUR_LOWER ) {
			Stg_KSPSetOperators( s->ksp_1, s->schur, s->schur, SAME_NONZERO_PATTERN );
			Stg_KSPSetOperators( s->ksp_2, A22, A22, SAME_NONZERO_PATTERN );
		}
		
	}
	
	/* do this once */
	if( pc->setfromoptionscalled ) {
		KSPSetFromOptions( s->ksp_1 );
		KSPSetFromOptions( s->ksp_2 );
	}
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode PCApply_BSchur_UPPER(PC pc,Vec x,Vec y)
{
	PC_BSchur s = (PC_BSchur)pc->data;
	Mat A11,A12,A21,A22;
	Vec x1,x2,y1,y2;
	Mat B;
	Vec t1;
	
	
	
	PetscFunctionBegin;
	
	Stg_PCGetOperators( pc, 0, &B, 0 );
	MatBlockGetSubMatrix( B, 0,0, &A11 );
	MatBlockGetSubMatrix( B, 0,1, &A12 );
	MatBlockGetSubMatrix( B, 1,0, &A21 );
	MatBlockGetSubMatrix( B, 1,1, &A22 );
	
	VecBlockGetSubVector( x, 0, &x1 );
	VecBlockGetSubVector( x, 1, &x2 );
	VecBlockGetSubVector( y, 0, &y1 );
	VecBlockGetSubVector( y, 1, &y2 );
	
	VecDuplicate( x1, &t1 );
	
	
	if( s->x_prime == PETSC_NULL) {
		VecDuplicate( x2, &s->x_prime );
	}
	
	/* form x prime : A21 inv(A11) x1 - x2 */
	if( s->block_factorisation == PETSC_TRUE ) {
		KSPSolve( s->ksp_1, x1, t1 );    /* t1 ~ inv(A11) x1 */
		MatMult( A21, t1, s->x_prime );  /* x' = A21 t1 */
		VecAXPY( s->x_prime, -1.0, x2 ); /* x' = x' - x2 */
	}
	else {
		VecCopy( x2, s->x_prime );
	}
	
	
	/* (2,2) solve */
	KSPSolve( s->ksp_2, s->x_prime, y2 );    /* y2 ~ inv(S) x' */
	
	/* form rhs */
	MatMult( A12, y2, t1 );       /* t1 = A12 y2 */
	VecScale( t1, -1.0 );         /* t1 = -t1 */
	VecAXPY( t1, 1.0, x1 );       /* t1 = t1 + x1 */
	/* (1,1) solve */
	KSPSolve( s->ksp_1, t1, y1 ); /* y1 ~ inv(A11) t1 */
	
	
	Stg_VecDestroy( &t1);
	
	PetscFunctionReturn(0);
}
PetscErrorCode PCApplyTranspose_BSchur_UPPER(PC pc,Vec x,Vec y)
{
    //PC_BSchur s = (PC_BSchur)pc->data;
	
	PetscFunctionBegin; 
	
	PetscFunctionReturn(0);
}

PetscErrorCode PCApply_BSchur_LOWER(PC pc,Vec x,Vec y)
{
	PC_BSchur s = (PC_BSchur)pc->data;
	Mat A11,A12,A21,A22;
	Vec x1,x2,y1,y2;
	Mat B;
	Vec t2;
	
	
	
	PetscFunctionBegin;
	
	Stg_PCGetOperators( pc, 0, &B, 0 );
	MatBlockGetSubMatrix( B, 0,0, &A11 );
	MatBlockGetSubMatrix( B, 0,1, &A12 );
	MatBlockGetSubMatrix( B, 1,0, &A21 );
	MatBlockGetSubMatrix( B, 1,1, &A22 );
	
	VecBlockGetSubVector( x, 0, &x1 );
	VecBlockGetSubVector( x, 1, &x2 );
	VecBlockGetSubVector( y, 0, &y1 );
	VecBlockGetSubVector( y, 1, &y2 );
	
	VecDuplicate( x2, &t2 );
	
	
	if( s->x_prime == PETSC_NULL) {
		VecDuplicate( x1, &s->x_prime );
	}
	
	/* form x prime : A12 inv(A22) x2 - x1 */
	if( s->block_factorisation == PETSC_TRUE ) {
		KSPSolve( s->ksp_2, x2, t2 );    /* t2 ~ inv(A22) x2 */
		MatMult( A12, t2, s->x_prime );  /* x' = A12 t2 */
		VecAXPY( s->x_prime, -1.0, x1 ); /* x' = x' - x1 */
	}
	else {
		VecCopy( x1, s->x_prime );
	}
	
	
	/* (1,1) solve */
	KSPSolve( s->ksp_1, s->x_prime, y1 );    /* y1 ~ inv(S) x' */
	
	/* form rhs */
	MatMult( A21, y1, t2 );       /* t2 = A21 y1 */
	VecScale( t2, -1.0 );         /* t2 = -t2 */
	VecAXPY( t2, 1.0, x2 );       /* t2 = t2 + x2 */
	/* (2,2) solve */
	KSPSolve( s->ksp_2, t2, y2 ); /* y2 ~ inv(A22) t2 */
	
	
	Stg_VecDestroy( &t2);
	
	PetscFunctionReturn(0);
}




PetscErrorCode PCSetFromOptions_BSchur(PC pc)
{
	PC_BSchur      s = (PC_BSchur)pc->data;
	PetscTruth     flg,set;
	PCBSchurType   type;
    PetscErrorCode ierr;
	const char *pc_prefix;
	char *option_name;
	
	PetscFunctionBegin;
	
	
	PCGetOptionsPrefix( pc, &pc_prefix );
	
	
	
	
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD, PETSC_NULL, "PC BSchur Options", "PC");CHKERRQ(ierr);
	
	/* type */
	if( asprintf( &option_name, "-%sbschur_type", pc_prefix ) > 0 ){
      PetscOptionsEnum(option_name,"Specifies the block structure of the schur preconditioner","PCBSchurSetType",
                       PCBSchurTypes, (PetscEnum)s->application_type,
                       (PetscEnum*)&type, &flg );
      if (flg) {
		PCBSchurSetType( pc, type );
      }
      free(option_name);
	}else{
      PetscPrintf( PETSC_COMM_SELF, "  Failed to create prefix for BSchur PC");
    }
	
	/* explicit */
	if( asprintf( &option_name, "-%sbschur_explicit", pc_prefix ) > 0 ){
      PetscOptionsTruth(option_name,"Specifies that we will form an explicit schur complement using a diagonal approximation","PCBSchurSetExplicit",
                        PETSC_FALSE, &flg, &set );
      if (flg) {
		s->explicit_operator = PETSC_TRUE;
      }
      free(option_name);
	}else{
      PetscPrintf( PETSC_COMM_SELF, "  Failed to create option_name for explicit");
    }
	
	/* use true rhs */
	if( asprintf( &option_name, "-%sbschur_use_factored_rhs", pc_prefix ) > 0 ){
      PetscOptionsTruth(option_name,"Specifies that the rhs used will be consistent with that obtained from doing block factorisation","PCBSchurSetUseFactoredRhs",
                        PETSC_FALSE, &flg, &set );
      if (flg) {
		s->block_factorisation = PETSC_TRUE;
      }
      free(option_name);
	}else{
      PetscPrintf( PETSC_COMM_SELF, "  Failed to create option_name for bschur_use_factored_rhs");
    }
	
	
	ierr = PetscOptionsEnd();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



PetscErrorCode PCView_BSchur(PC pc,PetscViewer viewer)
{
	PC_BSchur       s = (PC_BSchur)pc->data;
	PetscTruth isascii;
		
	
	PetscFunctionBegin; 
	
	Stg_PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&isascii);
	if (isascii) {
	
		PetscViewerASCIIPushTab( viewer );
		PetscViewerASCIIPrintf(viewer,"BSchur: type %s \n", PCBSchurTypes[ s->application_type ] );
		
		if(s->explicit_operator==PETSC_TRUE) {
			PetscViewerASCIIPrintf(viewer,"BSchur:  using explicit operator for schur complement \n" );
		}
		else {
			PetscViewerASCIIPrintf(viewer,"BSchur:  using matrix-free operator for schur complement \n" );
		}
		
		if(s->block_factorisation==PETSC_TRUE) {
			PetscViewerASCIIPrintf(viewer,"BSchur:  using modified rhs vector consistent with block factorisation \n" );
		}
		else {
			PetscViewerASCIIPrintf(viewer,"BSchur:  using original rhs vector \n" );
		}
		
		PetscViewerASCIIPrintf(viewer,"BSchur:  KSP (1,1) \n" );
		PetscViewerASCIIPushTab( viewer );
		KSPView( s->ksp_1, viewer );
		PetscViewerASCIIPopTab( viewer );
		
		PetscViewerASCIIPrintf(viewer,"BSchur:  KSP (2,2) \n" );
		PetscViewerASCIIPushTab( viewer );
		KSPView( s->ksp_2, viewer );
		PetscViewerASCIIPopTab( viewer );
		
		
		PetscViewerASCIIPopTab( viewer );
	}
	PetscFunctionReturn(0);
}




