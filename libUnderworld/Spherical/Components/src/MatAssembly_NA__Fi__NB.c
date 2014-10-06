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
#include "MatAssembly_NA__Fi__NB.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type MatAssembly_NA__Fi__NB_Type = "MatAssembly_NA__Fi__NB";

/* Creation implementation / Virtual constructor */
MatAssembly_NA__Fi__NB* _MatAssembly_NA__Fi__NB_New(  MATRIXASSEMBLYTERM_NA__NB__F_DEFARGS  )
{
   MatAssembly_NA__Fi__NB* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(MatAssembly_NA__Fi__NB) );
   self = (MatAssembly_NA__Fi__NB*) _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_PASSARGS  );

   /* Virtual info */

   return self;
}

void _MatAssembly_NA__Fi__NB_Init( void* matrixTerm, PpcManager* ppcManager, int grad_rho, int rho ) {
   MatAssembly_NA__Fi__NB* self = (MatAssembly_NA__Fi__NB*)matrixTerm;

   self->errorStream   = Journal_Register( Error_Type, (Name)self->name  );
   self->ppcManager    = ppcManager;
   self->grad_rho = grad_rho;
   self->rho = rho;

}

void _MatAssembly_NA__Fi__NB_Delete( void* matrixTerm ) {
   MatAssembly_NA__Fi__NB* self = (MatAssembly_NA__Fi__NB*)matrixTerm;

   _StiffnessMatrixTerm_Delete( self );
}

void _MatAssembly_NA__Fi__NB_Print( void* matrixTerm, Stream* stream ) {
   MatAssembly_NA__Fi__NB* self = (MatAssembly_NA__Fi__NB*)matrixTerm;

   _StiffnessMatrixTerm_Print( self, stream );

   /* General info */
}

void* _MatAssembly_NA__Fi__NB_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                                 _sizeOfSelf = sizeof(MatAssembly_NA__Fi__NB);
   Type                                                         type = MatAssembly_NA__Fi__NB_Type;
   Stg_Class_DeleteFunction*                                 _delete = _MatAssembly_NA__Fi__NB_Delete;
   Stg_Class_PrintFunction*                                   _print = _MatAssembly_NA__Fi__NB_Print;
   Stg_Class_CopyFunction*                                     _copy = NULL;
   Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _MatAssembly_NA__Fi__NB_DefaultNew;
   Stg_Component_ConstructFunction*                       _construct = _MatAssembly_NA__Fi__NB_AssignFromXML;
   Stg_Component_BuildFunction*                               _build = _MatAssembly_NA__Fi__NB_Build;
   Stg_Component_InitialiseFunction*                     _initialise = _MatAssembly_NA__Fi__NB_Initialise;
   Stg_Component_ExecuteFunction*                           _execute = _MatAssembly_NA__Fi__NB_Execute;
   Stg_Component_DestroyFunction*                           _destroy = _MatAssembly_NA__Fi__NB_Destroy;
   StiffnessMatrixTerm_AssembleElementFunction*     _assembleElement = _MatAssembly_NA__Fi__NB_AssembleElement;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*)_MatAssembly_NA__Fi__NB_New(  MATRIXASSEMBLYTERM_NA__NB__F_PASSARGS  );
}

void _MatAssembly_NA__Fi__NB_AssignFromXML( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
   MatAssembly_NA__Fi__NB*            self             = (MatAssembly_NA__Fi__NB*)matrixTerm;
   PpcManager*                ppcManager  = NULL;
   int                        grad_rho, rho;

   /* Construct Parent */
   _StiffnessMatrixTerm_AssignFromXML( self, cf, data );

   /* The PpcManager */
   ppcManager = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
   if( !ppcManager  )
      ppcManager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data  );

   grad_rho = PpcManager_GetPpcFromDict( ppcManager, cf, self->name, (Dictionary_Entry_Key)"grad_rho", "grad_rho" );
   rho = PpcManager_GetPpcFromDict( ppcManager, cf, self->name, (Dictionary_Entry_Key)"rho", "rho" );

   _MatAssembly_NA__Fi__NB_Init( self, ppcManager, grad_rho, rho );
}

void _MatAssembly_NA__Fi__NB_Build( void* matrixTerm, void* data ) {
   MatAssembly_NA__Fi__NB*             self             = (MatAssembly_NA__Fi__NB*)matrixTerm;

    Stg_Component_Build( self->ppcManager, data, False );
   _StiffnessMatrixTerm_Build( self, data );
}

void _MatAssembly_NA__Fi__NB_Initialise( void* matrixTerm, void* data ) {
   MatAssembly_NA__Fi__NB* self = (MatAssembly_NA__Fi__NB*)matrixTerm;

   Stg_Component_Initialise( self->ppcManager, data, False );
   _StiffnessMatrixTerm_Initialise( self, data );
}

void _MatAssembly_NA__Fi__NB_Execute( void* matrixTerm, void* data ) {
   _StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _MatAssembly_NA__Fi__NB_Destroy( void* matrixTerm, void* data ) {
   _StiffnessMatrixTerm_Destroy( matrixTerm, data );
}

void _MatAssembly_NA__Fi__NB_AssembleElement(
      void*                                              matrixTerm,
      StiffnessMatrix*                                   stiffnessMatrix,
      Element_LocalIndex                                 lElement_I,
      SystemLinearEquations*                             sle,
      FiniteElementContext*                              context,
      double**                                           elStiffMat )
{
   MatAssembly_NA__Fi__NB*   self = (MatAssembly_NA__Fi__NB*)matrixTerm;
   Swarm*                              swarm        = self->integrationSwarm;
   FeVariable*                         rowFeVar     = stiffnessMatrix->rowVariable;
   FeVariable*                         colFeVar     = stiffnessMatrix->columnVariable;
   Dimension_Index                     dim          = stiffnessMatrix->dim;
   int                                 rowDofs = rowFeVar->fieldComponentCount; // number of dofs per row node
   int                                 colDofs = colFeVar->fieldComponentCount; // number of dofs per row node
   IntegrationPoint*                   currIntegrationPoint;
   double*                             xi;
   double                              weight;
   Particle_InCellIndex                cParticle_I, cellParticleCount;
   Index                               rowNodes; // number of row nodes per element
   Index                               colNodes; // number of col nodes per element
   Index                               A,B;
   Index                               i;
   double                              detJac;
   double                              gradRho_rtp[3], gradRho_xyz[3], rho, xyz[3];
   Cell_Index                          cell_I;
   ElementType*                        rowElementType, *colElementType;
   double                              N[27], M[6];

   /* Set the element type */
   rowElementType = FeMesh_GetElementType( rowFeVar->feMesh, lElement_I ); rowNodes = rowElementType->nodeCount;
   colElementType = FeMesh_GetElementType( colFeVar->feMesh, lElement_I ); colNodes = colElementType->nodeCount;

   cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
   cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

   for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {

      currIntegrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );

      xi = currIntegrationPoint->xi;
      weight = currIntegrationPoint->weight;

      /* Calculate Determinant of Jacobian and Shape Functions */
      detJac = ElementType_JacobianDeterminant( colElementType, colFeVar->feMesh, lElement_I, xi, dim );
      ElementType_EvaluateShapeFunctionsAt( rowElementType, xi, M );
      ElementType_EvaluateShapeFunctionsAt( colElementType, xi, N ); 

      /* evaluate ppc function */
      PpcManager_Get( self->ppcManager, lElement_I, currIntegrationPoint, self->grad_rho, &(gradRho_rtp[0]) );
      PpcManager_Get( self->ppcManager, lElement_I, currIntegrationPoint, self->rho, &rho );
      FeMesh_CoordLocalToGlobal( colFeVar->feMesh, lElement_I, xi, xyz );
      /*
      gradRho_rtp[0]=gradRho_rtp[0]/rho; gradRho_rtp[1]=0; gradRho_rtp[2]=0;
      Spherical_VectorRTP2XYZ( gradRho_rtp, xyz, 2, gradRho_xyz );
      */
      gradRho_xyz[0]=0;
      gradRho_xyz[1]=gradRho_rtp[1]/rho;

      for( A=0; A<rowNodes; A++ ) {
         for( B=0; B<colNodes; B++ ) {
            for ( i = 0; i < colDofs ; i++ ) {
               elStiffMat[rowDofs*A][colDofs*B+i] += detJac * weight * M[A] * gradRho_xyz[i] * N[B] ;
            }
         }
      }
   }
}


