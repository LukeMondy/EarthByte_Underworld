/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	AuScope - http://www.auscope.org
**	Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"
#include <petscblaslapack.h>

/* Textual name of this class */
const Type SphericalStiffnessMatrix_Type = "SphericalStiffnessMatrix";

/** First part of name for build entry point */
static const char	SphericalStiffnessMatrix_assembleStiffnessMatrixStr[] = "assembleSphericalStiffnessMatrix";


void* SphericalStiffnessMatrix_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                                               _sizeOfSelf = sizeof(SphericalStiffnessMatrix);
   Type                                                                       type = SphericalStiffnessMatrix_Type;
   Stg_Class_DeleteFunction*                                               _delete = _StiffnessMatrix_Delete;
   Stg_Class_PrintFunction*                                                 _print = _StiffnessMatrix_Print;
   Stg_Class_CopyFunction*                                                   _copy = _StiffnessMatrix_Copy;
   Stg_Component_DefaultConstructorFunction*                   _defaultConstructor = StiffnessMatrix_DefaultNew;
   Stg_Component_ConstructFunction*                                     _construct = _SphericalStiffnessMatrix_AssignFromXML;
   Stg_Component_BuildFunction*                                             _build = _SphericalStiffnessMatrix_Build;
   Stg_Component_InitialiseFunction*                                   _initialise = _StiffnessMatrix_Initialise;
   Stg_Component_ExecuteFunction*                                         _execute = _StiffnessMatrix_Execute;
   Stg_Component_DestroyFunction*                                         _destroy = _SphericalStiffnessMatrix_Destroy;
   Bool                                                                   initFlag = False;
   StiffnessMatrix_CalculateNonZeroEntriesFunction*       _calculateNonZeroEntries = StiffnessMatrix_CalcNonZeros;
   void*                                                               rowVariable = NULL;
   void*                                                            columnVariable = NULL;
   void*                                                                       rhs = NULL;
   Stg_Component*                                               applicationDepInfo = NULL;
   Dimension_Index                                                             dim = 0;
   Bool                                                                isNonLinear = False;
   Bool                                              allowZeroElementContributions = False;
   void*                                                       entryPoint_Register = NULL;
   MPI_Comm                                                                   comm = 0;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return _SphericalStiffnessMatrix_New(  STIFFNESSMATRIX_PASSARGS  );
}


SphericalStiffnessMatrix* _SphericalStiffnessMatrix_New(  SPHERICALSTIFFNESSMATRIX_DEFARGS  )
{
   SphericalStiffnessMatrix*	self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(SphericalStiffnessMatrix) );
   /* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
   /* This means that any values of these parameters that are passed into this function are not passed onto the parent function
   and so should be set to ZERO in any children of this class. */
   nameAllocationType = NON_GLOBAL;

   self = (SphericalStiffnessMatrix*)_StiffnessMatrix_New(  STIFFNESSMATRIX_PASSARGS  );

   if( initFlag )   // 'initFlag' is confusing TODO: investigate it's removal
   {
      /* Only the EP self->assembleStiffnessMatrix is redefinied in this function */
      _SphericalStiffnessMatrix_Init( self );
   }

   return self;
}

void _SphericalStiffnessMatrix_Init(
   SphericalStiffnessMatrix*  self)
{

   if( self->assembleStiffnessMatrix )
   {
      // Override the behaviour setup by _StiffnessMatrix_Init
      EP_ReplaceAll( self->assembleStiffnessMatrix, _SphericalStiffnessMatrix_NewAssemble );
   }

   self->rotMat = NULL;
   self->tmpMat = NULL;
}

void _SphericalStiffnessMatrix_AssignFromXML( void* stiffnessMatrix, Stg_ComponentFactory* cf, void* data )
{
   SphericalStiffnessMatrix *self = (SphericalStiffnessMatrix*)stiffnessMatrix;

   _StiffnessMatrix_AssignFromXML( self, cf, data );

   _SphericalStiffnessMatrix_Init( self );
}

void _SphericalStiffnessMatrix_Build( void* stiffnessMatrix, void* data )
{
   SphericalStiffnessMatrix* self = (SphericalStiffnessMatrix*)stiffnessMatrix;

   _StiffnessMatrix_Build( stiffnessMatrix, data );

   // alloc some memory to these guys
   self->rotMat = ReallocArray( self->rotMat, double, (24*24) );
   self->tmpMat = ReallocArray( self->tmpMat, double, (24*24) );
}


void _SphericalStiffnessMatrix_Destroy( void* stiffnessMatrix, void* data )
{
   SphericalStiffnessMatrix* self = (SphericalStiffnessMatrix*)stiffnessMatrix;

   FreeArray( self->rotMat );
   FreeArray( self->tmpMat );

   _StiffnessMatrix_Destroy( stiffnessMatrix, data );


}

void _SphericalStiffnessMatrix_NewAssemble( void* stiffnessMatrix, Bool removeBCs, void* _sle, void* _context )
{
   static const double one = 1.0;
   SphericalStiffnessMatrix*		self = (SphericalStiffnessMatrix*)stiffnessMatrix;
   SystemLinearEquations*		sle = (SystemLinearEquations*)_sle;
   FeVariable			*rowVar, *colVar;
   FeMesh				*rowMesh, *colMesh;
   FeEquationNumber		*rowEqNum, *colEqNum;
   DofLayout			*rowDofs, *colDofs;
   unsigned			nRowEls;
   int			nRowNodes, *rowNodes;
   int			nColNodes, *colNodes;
   unsigned			maxDofs, maxRCDofs, nDofs, nRowDofs, nColDofs;
   double**			elStiffMat;
   double*				bcVals;
   Mat                             matrix = self->matrix;
   Vec				vector, transVector;
   int nRowNodeDofs, nColNodeDofs;
   int rowInd, colInd;
   double bc;
   unsigned			e_i, n_i, dof_i, n_j, dof_j;

   assert( self && Stg_CheckType( self, StiffnessMatrix ) );

   rowVar = self->rowVariable;
   colVar = self->columnVariable ? self->columnVariable : rowVar;
   rowEqNum = rowVar->eqNum;
   colEqNum = colVar->eqNum;
   rowMesh = rowVar->feMesh;
   colMesh = colVar->feMesh;
   rowDofs = rowVar->dofLayout;
   colDofs = colVar->dofLayout;
   nRowEls = FeMesh_GetElementLocalSize( rowMesh );
   assert( (rowVar == colVar) ? !self->transRHS : 1 );

   //matrix = self->matrix;
   vector = self->rhs ? self->rhs->vector : NULL;
   transVector = self->transRHS ? self->transRHS->vector : NULL;
   elStiffMat = NULL;
   bcVals = NULL;
   maxDofs = 0;

   /* Begin assembling each element. */
   for( e_i = 0; e_i < nRowEls; e_i++ )
   {
      FeMesh_GetElementNodes( rowMesh, e_i, self->rowInc );
      nRowNodes = IArray_GetSize( self->rowInc );
      rowNodes = IArray_GetPtr( self->rowInc );
      FeMesh_GetElementNodes( colMesh, e_i, self->colInc );
      nColNodes = IArray_GetSize( self->colInc );
      colNodes = IArray_GetPtr( self->colInc );

      /* Do we need more space to assemble this element? */
      nRowDofs = 0;
      for( n_i = 0; n_i < nRowNodes; n_i++ )
         nRowDofs += rowDofs->dofCounts[rowNodes[n_i]];
      nColDofs = 0;
      for( n_i = 0; n_i < nColNodes; n_i++ )
         nColDofs += colDofs->dofCounts[colNodes[n_i]];
      nDofs = nRowDofs * nColDofs;
      self->nRowDofs = nRowDofs;

      self->nColDofs = nColDofs;
      if( nDofs > maxDofs )
      {
         maxRCDofs = (nRowDofs > nColDofs) ? nRowDofs : nColDofs;
         elStiffMat = ReallocArray2D( elStiffMat, double, nRowDofs, nColDofs );
         bcVals = ReallocArray( bcVals, double, maxRCDofs );
         /* allocate the memory for a rotation matrix, only used if FeVariables involved
         	 have non axis aligned boundary conditions (nonAABCs) set to True */
         /*
         if( rowVar->nonAABCs ) {
         	self->rotMat = realloc( self->rotMat, maxRCDofs*maxRCDofs*sizeof(double) );
         	self->tmpMat = realloc( self->tmpMat, maxRCDofs*maxRCDofs*sizeof(double) );
         }  else if( colVar->nonAABCs ) {
         	self->rotMat = realloc( self->rotMat, maxRCDofs*maxRCDofs*sizeof(double) );
         	self->tmpMat = realloc( self->tmpMat, maxRCDofs*maxRCDofs*sizeof(double) );
         }
         */
         maxDofs = nDofs;
         self->elStiffMat = elStiffMat;
         self->bcVals = bcVals;
      }

      /* Assemble the element. */
      memset( elStiffMat[0], 0, nDofs * sizeof(double) );
      SphericalStiffnessMatrix_AssembleElement( self, e_i, sle, _context, elStiffMat );

      /* Correct for BCs providing I'm not keeping them in. */
      if( vector )
      {
         memset( bcVals, 0, nRowDofs * sizeof(double) );

         rowInd = 0;
         for( n_i = 0; n_i < nRowNodes; n_i++ )
         {
            nRowNodeDofs = rowDofs->dofCounts[rowNodes[n_i]];
            for( dof_i = 0; dof_i < nRowNodeDofs; dof_i++ )
            {
               if( !FeVariable_IsBC( rowVar, rowNodes[n_i], dof_i ) )
               {
                  colInd = 0;
                  for( n_j = 0; n_j < nColNodes; n_j++ )
                  {
                     nColNodeDofs = colDofs->dofCounts[colNodes[n_j]];
                     for( dof_j = 0; dof_j < nColNodeDofs; dof_j++ )
                     {
                        if( FeVariable_IsBC( colVar, colNodes[n_j], dof_j ) )
                        {
                           bc = DofLayout_GetValueDouble( colDofs, colNodes[n_j], dof_j );
                           bcVals[rowInd] -= bc * elStiffMat[rowInd][colInd];
                        }
                        colInd++;
                     }
                  }
               }
               rowInd++;
            }
         }

         VecSetValues( vector, nRowDofs, (int*)rowEqNum->locationMatrix[e_i][0], bcVals, ADD_VALUES );
      }
      if( transVector )
      {
         memset( bcVals, 0, nColDofs * sizeof(double) );

         colInd = 0;
         for( n_i = 0; n_i < nColNodes; n_i++ )
         {
            nColNodeDofs = colDofs->dofCounts[colNodes[n_i]];
            for( dof_i = 0; dof_i < nColNodeDofs; dof_i++ )
            {
               if( !FeVariable_IsBC( colVar, colNodes[n_i], dof_i ) )
               {
                  rowInd = 0;
                  for( n_j = 0; n_j < nRowNodes; n_j++ )
                  {
                     nRowNodeDofs = rowDofs->dofCounts[rowNodes[n_j]];
                     for( dof_j = 0; dof_j < nRowNodeDofs; dof_j++ )
                     {
                        if( FeVariable_IsBC( rowVar, rowNodes[n_j], dof_j ) )
                        {
                           bc = DofLayout_GetValueDouble( rowDofs, rowNodes[n_j], dof_j );
                           bcVals[colInd] -= bc * elStiffMat[rowInd][colInd];
                        }
                        rowInd++;
                     }
                  }
               }
               colInd++;
            }
         }

         VecSetValues( transVector, nColDofs, (int*)colEqNum->locationMatrix[e_i][0], bcVals, ADD_VALUES );
      }

      /* If keeping BCs in, zero corresponding entries in the element stiffness matrix. */
      if( !rowEqNum->removeBCs || !colEqNum->removeBCs )
      {
         rowInd = 0;
         for( n_i = 0; n_i < nRowNodes; n_i++ )
         {
            nRowNodeDofs = rowDofs->dofCounts[rowNodes[n_i]];
            for( dof_i = 0; dof_i < nRowNodeDofs; dof_i++ )
            {
               if( FeVariable_IsBC( rowVar, rowNodes[n_i], dof_i ) )
               {
                  memset( elStiffMat[rowInd], 0, nColDofs * sizeof(double) );
               }
               else
               {
                  colInd = 0;
                  for( n_j = 0; n_j < nColNodes; n_j++ )
                  {
                     nColNodeDofs = colDofs->dofCounts[colNodes[n_j]];
                     for( dof_j = 0; dof_j < nColNodeDofs; dof_j++ )
                     {
                        if( FeVariable_IsBC( colVar, colNodes[n_j], dof_j ) )
                           elStiffMat[rowInd][colInd] = 0.0;
                        colInd++;
                     }
                  }
               }
               rowInd++;
            }
         }
      }

      /* Add to stiffness matrix. */
      MatSetValues( matrix,
                    nRowDofs, (int*)rowEqNum->locationMatrix[e_i][0],
                    nColDofs, (int*)colEqNum->locationMatrix[e_i][0],
                    elStiffMat[0], ADD_VALUES );
   }

   FreeArray( elStiffMat );
   FreeArray( bcVals );

   /* If keeping BCs in and rows and columnns use the same variable, put ones in all BC'd diagonals. */
   if( !colEqNum->removeBCs && rowVar == colVar )
   {
      for( n_i = 0; n_i < FeMesh_GetNodeLocalSize( colMesh ); n_i++ )
      {
         nColNodeDofs = colDofs->dofCounts[n_i];
         for( dof_i = 0; dof_i < nColNodeDofs; dof_i++ )
         {
            if( FeVariable_IsBC( colVar, n_i, dof_i ) )
            {
               MatSetValues( self->matrix,
                             1, colEqNum->destinationArray[n_i] + dof_i,
                             1, colEqNum->destinationArray[n_i] + dof_i,
                             (double*)&one, ADD_VALUES );
            }
         }
      }
   }

   /* Reassemble the matrix and vectors. */
   MatAssemblyBegin( matrix, MAT_FINAL_ASSEMBLY );
   MatAssemblyEnd( matrix, MAT_FINAL_ASSEMBLY );
   if( vector )
   {
      VecAssemblyBegin( vector );
      VecAssemblyEnd( vector );
   }
   if( transVector)
   {
      VecAssemblyBegin( transVector );
      VecAssemblyEnd( transVector );
   }

   MatAssemblyBegin( matrix, MAT_FINAL_ASSEMBLY );
   MatAssemblyEnd( matrix, MAT_FINAL_ASSEMBLY );
}

void blasMatrixMult( double A[], double B[], int rowA, int colB, int colA, double *C )
{
   /*@
   	Performs [C] = [A]*[B],

   	where [A], [B], [C] are single arrays, with rowMajor ordering, so they can be vectors.
        [C] must be pre allocated
   @*/
   char n='N';
   PetscScalar zero=0.0;
   PetscScalar alpha=1.0;

   /* BLASgemm expects fortan memory chunks, so the problem is formed using
   	 the matrix transpose formula [C]^T = [B]^T * [A]^T */

   BLASgemm_( &n, &n,
              &colB, &rowA, &colA,
              &alpha, B, &colB,
              A, &colA, &zero,
              C, &colB );
}


void SphericalStiffnessMatrix_AssembleElement(
   void* stiffnessMatrix,
   Element_LocalIndex element_lI,
   SystemLinearEquations* sle,
   FiniteElementContext* context,
   double** elStiffMatToAdd )
{
   SphericalStiffnessMatrix* self = (SphericalStiffnessMatrix*) stiffnessMatrix;
   Index                   stiffnessMatrixTermCount  = Stg_ObjectList_Count( self->stiffnessMatrixTermList );
   Index                   stiffnessMatrixTerm_I;
   StiffnessMatrixTerm*    stiffnessMatrixTerm;

   FeVariable* rowVar=self->rowVariable;
   FeVariable* colVar=self->columnVariable;

   double* R = self->rotMat;
   double* tmp = self->tmpMat;
   Bool Reval=False; /* state flag if R has already been evaluated */

   for ( stiffnessMatrixTerm_I = 0 ; stiffnessMatrixTerm_I < stiffnessMatrixTermCount ; stiffnessMatrixTerm_I++ )
   {
      stiffnessMatrixTerm = (StiffnessMatrixTerm*) Stg_ObjectList_At( self->stiffnessMatrixTermList, stiffnessMatrixTerm_I );
      StiffnessMatrixTerm_AssembleElement( stiffnessMatrixTerm, (StiffnessMatrix*)self, element_lI, sle, context, elStiffMatToAdd );
   }

   /* check if we need to PRE MULTIPLY with rotation matrix, [R] -
   	 only if rowVariable has non axis aligned BCs and the element is on the mesh boundary */
   if( rowVar->nonAABCs && IndexSet_IsMember( rowVar->feMesh->bndElementSet, element_lI ) )
   {

      /*	Perform [tmp] = [R]^T * [elStiffMatToAdd]
      	 	but do it with BLAS (fortran column major ordered) memory layout
       		therefore compute: [tmp]^T = [elStiffMatToAdd]^T * [[R]^T]^T
       */

      int rowA = (self->nRowDofs > self->nColDofs) ? self->nRowDofs : self->nColDofs; // rows in [R]
      int colA = rowA;              // cols in [R]
      int colB = self->nColDofs;    // cols in [elStiffMatToAdd]

      PetscScalar one=1.0;
      PetscScalar zero=0.0;
      char n='N';
      char t='T';

      // initialise [R] and [tmp] memory
      memset(R,0,rowA*rowA*sizeof(double));
      memset(tmp,0,rowA*rowA*sizeof(double));

      // evaluate [R] for element_lI
      Reval = (self->dim == 3 ) ? 
	 SphericalStiffnessMatrix_EvaluateRotationMatrix3D( rowVar, element_lI, R ) : 
	 SphericalStiffnessMatrix_EvaluateRotationMatrix2D( rowVar, element_lI, R );

      // [tmp]^T = [elStiffMatToAdd]^T * [R]
      // only using BLASgemm_ [not blasMatrixMult()] because transposed matrices are used
      BLASgemm_( &n, &t, &colB, &rowA, &colA, &one,elStiffMatToAdd[0], &colB, R, &colA, &zero, tmp, &colB );

      // copy result into returned memory segment
      memcpy( elStiffMatToAdd[0], tmp, rowA*colB*sizeof(double) );
   }


   /* check if we need to POST MULTIPLY with rotation matrix -
   	 only if columnVariable has non axis aligned BCs and the element is on the mesh boundary */
   if( colVar->nonAABCs && IndexSet_IsMember( colVar->feMesh->bndElementSet, element_lI ) )
   {
      /*
         Perform [tmp] = [elStiffMatToAdd] * [R]
         but do it with BLAS (fortran column major ordered) memory layout
         therefore compute: [tmp]^T = [R]^T * [elStiffMatToAdd]^T 
       */

      int rowA = self->nRowDofs; // rows in [elStiffMatToAdd]
      int colB = (self->nRowDofs > self->nColDofs) ? self->nRowDofs : self->nColDofs; // colB in [R]
      int colA = self->nColDofs;   // cols in [R]

      // check if R has already been evaluated
      if( !Reval )
      {
         memset(R,0,colB*colB*sizeof(double));
         Reval = ( self->dim==3 ) ?
                 SphericalStiffnessMatrix_EvaluateRotationMatrix3D( colVar, element_lI, R ) :
                 SphericalStiffnessMatrix_EvaluateRotationMatrix2D( colVar, element_lI, R );
      }
      memset(tmp,0,colB*colB*sizeof(double));

      // tmp = [elStiffMatToAdd] * [R]
      blasMatrixMult( elStiffMatToAdd[0], R, rowA, colB, colA, tmp );

      memcpy( elStiffMatToAdd[0], tmp, rowA*colB*sizeof(double) );
   }
}

Bool SphericalStiffnessMatrix_EvaluateRotationMatrix2D( FeVariable* var, unsigned element_lI, double *rMat )
{
   /*@
   	calculate the rotation matrix for this element, assume it is already on the boundary of the domain

   TODO: not sure if nNdoes (for specifying how many nodes per element is really needed...

   @*/

   FeMesh *mesh = var->feMesh;
   IArray *inc = IArray_New();
   int *nodes, nNodes, node_I, start, offset, dof;
   double rot[9];
   /* Get nodes in element_lI */
   FeMesh_GetElementNodes( mesh, element_lI, inc );
   nNodes = IArray_GetSize( inc );
   nodes = IArray_GetPtr( inc );

   dof = var->fieldComponentCount;

   for( node_I=0; node_I<nNodes; node_I++ )
   {

      // start & offset numbers temp. place holder for building the elemental rotation matrix
      start = dof * node_I * ( dof*nNodes ) + dof*node_I;
      offset = dof*nNodes;

      if( !IndexSet_IsMember( mesh->bndNodeSet, nodes[node_I] ) )
      {
         /** if not a boundary node, build the Identity matrix for rotations, i.e no rotation */

         rMat[ start              ] = 1;
         rMat[ start +          1 ] = 0;

         rMat[ start + offset     ] = 0;
         rMat[ start + offset + 1 ] = 1;

      }
      else
      {
         /** if boundary node, build the proper rotation matrix */

         /** to calculate the rotation matrix */
         Spherical_Get_RotationMatrixIJK( (Mesh*)mesh, nodes[node_I], rot );

         rMat[ start              ] = rot[0];
         rMat[ start +          1 ] = rot[1];

         rMat[ start + offset     ] = rot[2];
         rMat[ start + offset + 1 ] = rot[3];
      }
   }
   Stg_Class_Delete( inc );
   return True;
}

Bool SphericalStiffnessMatrix_EvaluateRotationMatrix3D( FeVariable* var, unsigned element_lI, double *rMat )
{
   /*@
   	calculate the rotation matrix for this element, assume it is already on the boundary of the domain

   TODO: not sure if nNdoes (for specifying how many nodes per element is really needed...

   @*/

   FeMesh *mesh = var->feMesh;
   IArray *inc = IArray_New();
   int nNodes = FeMesh_GetElementNodeSize(mesh, element_lI);
   int *nodes, node_I, start, offset, dof;
   double rot[9];
   /* Get nodes in element_lI */
   FeMesh_GetElementNodes( mesh, element_lI, inc );
   nNodes = IArray_GetSize( inc );
   nodes = IArray_GetPtr( inc );

   dof = var->fieldComponentCount;

   for( node_I=0; node_I<nNodes; node_I++ )
   {

      // start & offset numbers temp. place holder for building the elemental rotation matrix
      start = dof * node_I * ( dof*nNodes ) + dof*node_I;
      offset = dof*nNodes;

      if( !IndexSet_IsMember( mesh->bndNodeSet, nodes[node_I] ) )
      {
         /** if not a boundary node, build the Identity matrix for rotations, i.e no rotation */

         rMat[ start                ] = 1;
         rMat[ start +            1 ] = 0;
         rMat[ start +            2 ] = 0;

         rMat[ start +   offset     ] = 0;
         rMat[ start +   offset + 1 ] = 1;
         rMat[ start +   offset + 2 ] = 0;

         rMat[ start + 2*offset     ] = 0;
         rMat[ start + 2*offset + 1 ] = 0;
         rMat[ start + 2*offset + 2 ] = 1;
      }
      else
      {
         /** if boundary node, build the proper rotation matrix */

         /** to calculate the rotation matrix */
         Spherical_Get_RotationMatrixIJK( (Mesh*)mesh, nodes[node_I], rot );

         rMat[ start                ] = rot[0];
         rMat[ start +            1 ] = rot[1];
         rMat[ start +            2 ] = rot[2];

         rMat[ start + 1*offset     ] = rot[3];
         rMat[ start + 1*offset + 1 ] = rot[4];
         rMat[ start + 1*offset + 2 ] = rot[5];

         rMat[ start + 2*offset     ] = rot[6];
         rMat[ start + 2*offset + 1 ] = rot[7];
         rMat[ start + 2*offset + 2 ] = rot[8];
      }
   }
   Stg_Class_Delete( inc );
   return True;
}

