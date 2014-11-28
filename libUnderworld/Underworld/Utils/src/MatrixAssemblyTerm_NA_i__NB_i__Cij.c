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


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "MatrixAssemblyTerm_NA_i__NB_i__Cij.h"

#include<petsc.h>
#include<petscblaslapack.h>
#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type MatrixAssemblyTerm_NA_i__NB_i__Cij_Type = "MatrixAssemblyTerm_NA_i__NB_i__Cij";

/* Creation implementation / Virtual constructor */
MatrixAssemblyTerm_NA_i__NB_i__Cij* _MatrixAssemblyTerm_NA_i__NB_i__Cij_New(  MATRIXASSEMBLYTERM_NA_I__NB_I__F_DEFARGS  )
{
   MatrixAssemblyTerm_NA_i__NB_i__Cij* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(MatrixAssemblyTerm_NA_i__NB_i__Cij) );
   self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*) _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_PASSARGS  );

   /* Virtual info */

   return self;
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Init( void* matrixTerm, PpcManager* ppcManager, int functionLabel ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij* self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;

   self->errorStream = Journal_Register( Error_Type, (Name)self->name  );
   self->ppcManager    = ppcManager;
   self->functionLabel = functionLabel;
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Delete( void* matrixTerm ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij* self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;

   _StiffnessMatrixTerm_Delete( self );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Print( void* matrixTerm, Stream* stream ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij* self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;

   _StiffnessMatrixTerm_Print( self, stream );

   /* General info */
}

void* _MatrixAssemblyTerm_NA_i__NB_i__Cij_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                 _sizeOfSelf = sizeof(MatrixAssemblyTerm_NA_i__NB_i__Cij);
   Type                                                         type = MatrixAssemblyTerm_NA_i__NB_i__Cij_Type;
   Stg_Class_DeleteFunction*                                 _delete = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Delete;
   Stg_Class_PrintFunction*                                   _print = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Print;
   Stg_Class_CopyFunction*                                     _copy = NULL;
   Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _MatrixAssemblyTerm_NA_i__NB_i__Cij_DefaultNew;
   Stg_Component_ConstructFunction*                       _construct = _MatrixAssemblyTerm_NA_i__NB_i__Cij_AssignFromXML;
   Stg_Component_BuildFunction*                               _build = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Build;
   Stg_Component_InitialiseFunction*                     _initialise = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Initialise;
   Stg_Component_ExecuteFunction*                           _execute = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Execute;
   Stg_Component_DestroyFunction*                           _destroy = _MatrixAssemblyTerm_NA_i__NB_i__Cij_Destroy;
   StiffnessMatrixTerm_AssembleElementFunction*     _assembleElement = _MatrixAssemblyTerm_NA_i__NB_i__Cij_AssembleElement;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*)_MatrixAssemblyTerm_NA_i__NB_i__Cij_New(  MATRIXASSEMBLYTERM_NA_I__NB_I__F_PASSARGS  );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_AssignFromXML( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij*            self             = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;
   PpcManager*                ppcManager  = NULL;
   int                        functionLabel;

   /* Construct Parent */
   _StiffnessMatrixTerm_AssignFromXML( self, cf, data );

   /* The PpcManager */
   ppcManager = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
   if( !ppcManager  )
      ppcManager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data  );

   functionLabel = PpcManager_GetPpcFromDict( ppcManager, cf, self->name, (Dictionary_Entry_Key)"functionLabel", "functionLabel" );
   self->nonlinearTag = PpcManager_GetPpcFromDict( ppcManager, cf, self->name, (Dictionary_Entry_Key)"nonlinear", "" );

   _MatrixAssemblyTerm_NA_i__NB_i__Cij_Init( self, ppcManager, functionLabel );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Build( void* matrixTerm, void* data ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij*             self             = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;

    Stg_Component_Build( self->ppcManager, data, False );
   _StiffnessMatrixTerm_Build( self, data );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Initialise( void* matrixTerm, void* data ) {
   MatrixAssemblyTerm_NA_i__NB_i__Cij* self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;
   double nonlinear;
   int p_i;
   IntegrationPoint *particle;

   Stg_Component_Initialise( self->ppcManager, data, False );
   _StiffnessMatrixTerm_Initialise( self, data );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Execute( void* matrixTerm, void* data ) {
   _StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_Destroy( void* matrixTerm, void* data ) {
   _StiffnessMatrixTerm_Destroy( matrixTerm, data );
}

void blasMatrixMultAB( double *A, double *B, int rowA, int colA, int colB, double *C )
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

void blasMatrixMultAtB( double *A, double *B, int rowA, int colA, int colB, double *C )
{
   /*@
      Performs [C] = [A]^T*[B],

      where [A], [B], [C] are single arrays, with rowMajor ordering, so they can be vectors.
      [C] must be pre allocated
   @*/
   char n='N';
   char t='T';
   PetscScalar zero=0.0;
   PetscScalar alpha=1.0;

   /* BLASgemm expects fortan memory chunks, so the problem is formed using
   the matrix transpose formula [C]^T = [B]^T * [A] */

   BLASgemm_( &n, &t,
      &colB, &rowA, &colA,
      &alpha, B, &colB,
      A, &rowA, &zero,
      C, &colB );
}

void Assemble_DB( double eta, double d_dx, double d_dy, double d_dz, int dim, double* B ){
   /*@
      Assembles Strain rate operator.
   @*/

   assert(dim>0);
   assert(B);

   if( dim == 2 ) {
   /*
   [B] = [ d/dx,     0  ]
         [    0,  d/dy  ]
         [ d/dy,  d/dx  ]  */

      B[0]=2*eta*d_dx;  B[1]=0;
      B[2]=0;     B[3]=2*eta*d_dy;
      B[4]=eta*d_dy;  B[5]=eta*d_dx;
   } else {
   /*
   [B] = [ d/dx,     0,      0  ]
         [    0,  d/dy,      0  ]
         [    0,     0,   d/dx  ]
         [ d/dy,  d/dx,      0  ]
         [ d/dz,     0,   d/dx  ]
         [    0,  d/dz,   d/dy  ] */

      B[0]=d_dx;  B[1]=0;     B[2]=0;
      B[3]=0;     B[4]=d_dy;  B[5]=0;
      B[6]=0;     B[7]=0;     B[8]=d_dz;
      B[9]=d_dy;  B[10]=d_dx; B[11]=0;
      B[12]=d_dz; B[13]=0;    B[14]=d_dx;
      B[15]=0;    B[16]=d_dz; B[17]=d_dy;
   }
 
} 
void Assemble_B( double d_dx, double d_dy, double d_dz, int dim, double* B ){
   /*@
      Assembles Strain rate operator.
   @*/

   assert(dim>0);
   assert(B);

   if( dim == 2 ) {
   /*
   [B] = [ d/dx,     0  ]
         [    0,  d/dy  ]
         [ d/dy,  d/dx  ]  */

      B[0]=d_dx;  B[1]=0;
      B[2]=0;     B[3]=d_dy;
      B[4]=d_dy;  B[5]=d_dx;
   } else {
   /*
   [B] = [ d/dx,     0,      0  ]
         [    0,  d/dy,      0  ]
         [    0,     0,   d/dx  ]
         [ d/dy,  d/dx,      0  ]
         [ d/dz,     0,   d/dx  ]
         [    0,  d/dz,   d/dy  ] */

      B[0]=d_dx;  B[1]=0;     B[2]=0;
      B[3]=0;     B[4]=d_dy;  B[5]=0;
      B[6]=0;     B[7]=0;     B[8]=d_dz;
      B[9]=d_dy;  B[10]=d_dx; B[11]=0;
      B[12]=d_dz; B[13]=0;    B[14]=d_dx;
      B[15]=0;    B[16]=d_dz; B[17]=d_dy;
   }
 
}
#if 0
void _MatrixAssemblyTerm_NA_i__NB_i__Cij_AssembleElement(
      void*                                              matrixTerm,
      StiffnessMatrix*                                   stiffnessMatrix,
      Element_LocalIndex                                 lElement_I,
      SystemLinearEquations*                             sle,
      FiniteElementContext*                              context,
      double**                                           elStiffMat )
{

   /**
   ASSUMPTION: The row and column variables that define this matrix are the same
   **/

   MatrixAssemblyTerm_NA_i__NB_i__Cij*   self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;
   Swarm*                              swarm        = self->integrationSwarm;
   FeVariable*                         rowVar       = stiffnessMatrix->rowVariable;
   Dimension_Index                     dim          = stiffnessMatrix->dim;
   IntegrationPoint*                   currIntegrationPoint;
   double*                             xi;
   double                              weight;
   Particle_InCellIndex                cParticle_I, cellParticleCount;
   Index                               nodesPerEl;
   Index                               a_i,b_i;
   Index                               i;
   int                                 rowDof_i, colDof_i, tensorDofCount; 
   double                              d_dx, d_dy, d_dz;
   double**                            GNx;
   double                              detJac;
   double                              F;
   Cell_Index                          cell_I;
   ElementType*                        elementType;
   unsigned                            dofCount;
   double                              B_a[18], B_b[18];
   double                              D[64]; // see louis paper for definition, is the consitutive matrix
   double                              DB[18]; // see louis paper for definition, is the consitutive matrix
   double                              BtDB[9]; // see louis paper for definition, is the consitutive matrix

   // sanity check and initialise all small matrices
   assert( stiffnessMatrix->rowVariable == stiffnessMatrix->columnVariable );
   memset( B_a, 0 ,sizeof(double)*18);
   memset( B_b, 0 ,sizeof(double)*18);
   memset( D, 0 ,sizeof(double)*64);
   memset( DB, 0 ,sizeof(double)*18);
   memset( BtDB, 0 ,sizeof(double)*9);

   /* Set the element type */
   elementType = FeMesh_GetElementType( rowVar->feMesh, lElement_I );
   nodesPerEl = elementType->nodeCount;
   dofCount = rowVar->fieldComponentCount;

   if( nodesPerEl > self->max_nElNodes ) {
      /* reallocate */
      // always assume 3 component (x,y,z) make programming easier
      self->GNx = (double*)ReallocArray2D( self->GNx, double, 3, nodesPerEl );
      self->max_nElNodes = nodesPerEl;
   }
   GNx = self->GNx;

   cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
   cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

   for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {

      currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

      xi = currIntegrationPoint->xi;
      weight = currIntegrationPoint->weight;

      ElementType_ShapeFunctionsGlobalDerivs(
         elementType,
         rowVar->feMesh, lElement_I,
         xi, dim, &detJac, GNx );

      /* evaluate ppc function */
      PpcManager_Get( self->ppcManager, lElement_I, currIntegrationPoint, self->functionLabel, &F );
      /* Get Constitutive matrix. 2D 3x3, 3D, 6x6 */
      if( dim == 2 ) {
         tensorDofCount=3;
         D[0] = D[4] = 2*F*detJac*weight; // components in plane
         D[8] = F*detJac*weight; // shear components 
      } else {
         tensorDofCount=6;
         D[0] = D[7] = D[14] = D[21] = D[28] = 2*F*detJac*weight; 
         D[21] = D[28] = D[35] = F*detJac*weight;
      }

      for( a_i=0; a_i<nodesPerEl; a_i++ ) {
         rowDof_i = a_i*dofCount;

         // assemble discrete gradient operator for row var
         d_dx = GNx[I_AXIS][a_i];
         d_dy = GNx[J_AXIS][a_i];
         d_dz = GNx[K_AXIS][a_i];

         Assemble_B( d_dx, d_dy, d_dz, dim, B_a );

         /* K = [D][B] */
         blasMatrixMultAB( D, B_a, tensorDofCount, tensorDofCount, dofCount, DB );

         for( b_i=0; b_i<nodesPerEl; b_i++ ) {
            colDof_i = b_i*dofCount;

            // assemble discrete gradient operator for col var
            d_dx = GNx[I_AXIS][b_i];
            d_dy = GNx[J_AXIS][b_i];
            d_dz = GNx[K_AXIS][b_i];
            Assemble_B( d_dx, d_dy, d_dz, dim, B_b );

            /* K = [B]^T[D][B] */
            blasMatrixMultAtB( B_b, DB, dofCount, tensorDofCount, dofCount, BtDB );

            if( dim == 2 ) {
               elStiffMat[colDof_i  ][rowDof_i  ] += BtDB[0];
               elStiffMat[colDof_i  ][rowDof_i+1] += BtDB[1];
               elStiffMat[colDof_i+1][rowDof_i  ] += BtDB[2];
               elStiffMat[colDof_i+1][rowDof_i+1] += BtDB[3];
            }
            if( dim == 3 ) {
               elStiffMat[colDof_i  ][rowDof_i  ] += BtDB[0];
               elStiffMat[colDof_i  ][rowDof_i+1] += BtDB[1];
               elStiffMat[colDof_i  ][rowDof_i+2] += BtDB[2];
               elStiffMat[colDof_i+1][rowDof_i  ] += BtDB[3];
               elStiffMat[colDof_i+1][rowDof_i+1] += BtDB[4];
               elStiffMat[colDof_i+1][rowDof_i+2] += BtDB[5];
               elStiffMat[colDof_i+2][rowDof_i  ] += BtDB[6];
               elStiffMat[colDof_i+2][rowDof_i+1] += BtDB[7];
               elStiffMat[colDof_i+2][rowDof_i+2] += BtDB[8];
            }
         }
      }
   }
}
#endif

void _MatrixAssemblyTerm_NA_i__NB_i__Cij_AssembleElement(
      void*                                              matrixTerm,
      StiffnessMatrix*                                   stiffnessMatrix,
      Element_LocalIndex                                 lElement_I,
      SystemLinearEquations*                             sle,
      FiniteElementContext*                              context,
      double**                                           elStiffMat )
{

   /**
   ASSUMPTION: The row and column variables that define this matrix are the same
   **/

   MatrixAssemblyTerm_NA_i__NB_i__Cij*   self = (MatrixAssemblyTerm_NA_i__NB_i__Cij*)matrixTerm;
   Swarm*                              swarm        = self->integrationSwarm;
   FeVariable*                         rowVar       = stiffnessMatrix->rowVariable;
   Dimension_Index                     dim          = stiffnessMatrix->dim;
   IntegrationPoint*                   currIntegrationPoint;
   double*                             xi;
   double                              weight;
   Particle_InCellIndex                cParticle_I, cellParticleCount;
   Index                               nodesPerEl;
   Index                               a_i,b_i;
   int                                 rowDof_i, colDof_i;
   double                              d_dx, d_dy, d_dz;
   double**                            GNx;
   double                              detJac;
   double                              F;
   Cell_Index                          cell_I;
   ElementType*                        elementType;
   unsigned                            dofCount;
   double                              B_a[18], B_b[18];
   double                              D[64]; // see louis paper for definition, is the consitutive matrix
   double                              DB[18]; // see louis paper for definition, is the consitutive matrix
   double                              BtDB[9]; // see louis paper for definition, is the consitutive matrix

   // sanity check and initialise all small matrices
   assert( stiffnessMatrix->rowVariable == stiffnessMatrix->columnVariable );
   memset( B_a, 0 ,sizeof(double)*18);
   memset( B_b, 0 ,sizeof(double)*18);
   memset( D, 0 ,sizeof(double)*64);
   memset( DB, 0 ,sizeof(double)*18);
   memset( BtDB, 0 ,sizeof(double)*9);

   /* Set the element type */
   elementType = FeMesh_GetElementType( rowVar->feMesh, lElement_I );
   nodesPerEl = elementType->nodeCount;
   dofCount = rowVar->fieldComponentCount;

   if( nodesPerEl > self->max_nElNodes ) {
      /* reallocate */
      // always assume 3 component (x,y,z) make programming easier
      self->GNx = (double*)ReallocArray2D( self->GNx, double, 3, nodesPerEl );
      self->max_nElNodes = nodesPerEl;
   }
   GNx = self->GNx;

   cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
   cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

   for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {

      currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

      xi = currIntegrationPoint->xi;
      weight = currIntegrationPoint->weight;

      ElementType_ShapeFunctionsGlobalDerivs(
         elementType,
         rowVar->feMesh, lElement_I,
         xi, dim, &detJac, GNx );

      
      /* evaluate ppc function */
      PpcManager_Get( self->ppcManager, lElement_I, currIntegrationPoint, self->functionLabel, &F );

      for( a_i=0; a_i<nodesPerEl; a_i++ ) {
         rowDof_i = a_i*dofCount;

         // evaluate discrete gradient operator for row var
         d_dx = GNx[I_AXIS][a_i];
         d_dy = GNx[J_AXIS][a_i];
         d_dz = GNx[K_AXIS][a_i];

         /* Build [D][B]. size [tensorDofCount,dofCount] */
         Assemble_DB( F*detJac*weight, d_dx, d_dy, d_dz, dim, DB );

         for( b_i=0; b_i<nodesPerEl; b_i++ ) {
            colDof_i = b_i*dofCount;

            // evaluate discrete gradient operator for col var
            d_dx = GNx[I_AXIS][b_i];
            d_dy = GNx[J_AXIS][b_i];
            d_dz = GNx[K_AXIS][b_i];

            /* K = [B]^T[D][B] */
            if( dim == 2 ) {
               elStiffMat[colDof_i  ][rowDof_i  ] += d_dx*DB[0]+d_dy*DB[4];
               elStiffMat[colDof_i  ][rowDof_i+1] += d_dy*DB[5];
               elStiffMat[colDof_i+1][rowDof_i  ] += d_dx*DB[4];
               elStiffMat[colDof_i+1][rowDof_i+1] += d_dy*DB[3]+d_dx*DB[5];
            }
            /* Need to update
            if( dim == 3 ) {
               elStiffMat[colDof_i  ][rowDof_i  ] += BtDB[0];
               elStiffMat[colDof_i  ][rowDof_i+1] += BtDB[1];
               elStiffMat[colDof_i  ][rowDof_i+2] += BtDB[2];
               elStiffMat[colDof_i+1][rowDof_i  ] += BtDB[3];
               elStiffMat[colDof_i+1][rowDof_i+1] += BtDB[4];
               elStiffMat[colDof_i+1][rowDof_i+2] += BtDB[5];
               elStiffMat[colDof_i+2][rowDof_i  ] += BtDB[6];
               elStiffMat[colDof_i+2][rowDof_i+1] += BtDB[7];
               elStiffMat[colDof_i+2][rowDof_i+2] += BtDB[8];
            }
            */
         }
      }
   }
}


