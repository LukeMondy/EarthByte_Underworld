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
#include "private/mat/petscext_matimpls.h"
#include <petscversion.h>
#include "private/compat/petsccompat.h"

EXTERN_C_BEGIN

EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatCreate_Block(Mat);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatCreate_SymTrans(Mat);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatCreate_Schur( Mat A );
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatCreate_RestrictScatter( Mat A );

EXTERN PetscLogStage __PETScExt_MatSchur_MatMult_Stage;
EXTERN PetscLogStage __PETScExt_MatSchur_MatMultTranspose_Stage;
EXTERN PetscLogEvent __PETScExt_MatSchur_MatMults_Event;
EXTERN PetscLogEvent __PETScExt_MatSchur_A11_inverse_Event;
EXTERN PetscLogEvent __PETScExt_MatSchur_A22_inverse_Event;


EXTERN_C_END


PetscErrorCode PetscExtMatRegisterAll( void )
{
	Stg_MatRegister(	"block",			"src/mat/impls/block",				"MatCreate_Block",				MatCreate_Block 		);
	Stg_MatRegister(	"symtrans",			"src/mat/impls/symtrans",			"MatCreate_SymTrans",			MatCreate_SymTrans 	);	
	Stg_MatRegister(	"schur",			"src/mat/impls/schur",				"MatCreate_Schur",				MatCreate_Schur 		);	
	Stg_MatRegister(	"restrictscatter",	"src/mat/impls/restrictscatter",	"MatCreate_RestrictScatter",	MatCreate_RestrictScatter	);	
	
	/* init logging for schur */
	PetscLogStageRegister( "PETScExt_MatSchur_MatMult", &__PETScExt_MatSchur_MatMult_Stage );
	PetscLogStageRegister( "PETScExt_MatSchur_MatMultTranspose", &__PETScExt_MatSchur_MatMultTranspose_Stage );
	
	PetscLogEventRegister( "PETScExt_MatSchur_MatMults", 0, &__PETScExt_MatSchur_MatMults_Event );
	PetscLogEventRegister( "PETScExt_MatSchur_A11_inverse", 0, &__PETScExt_MatSchur_A11_inverse_Event );
	PetscLogEventRegister( "PETScExt_MatSchur_A22_inverse", 0, &__PETScExt_MatSchur_A22_inverse_Event );
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode PetscExtMatRegisterDestroyAll( void )
{
	PetscFunctionReturn(0);
}


