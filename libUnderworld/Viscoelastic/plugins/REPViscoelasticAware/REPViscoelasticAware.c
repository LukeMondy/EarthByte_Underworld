/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005-2010, Monash University
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice,
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright
**			notice, this list of conditions and the following disclaimer in the
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors
**			may be used to endorse or promote products derived from this software
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%  Louis.Moresi - Louis.Moresi@monash.edu
*%
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* This plugin will cause ... */
#include <mpi.h>
#include <string.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Viscoelastic/Base/Base.h>

/* call these global vars a disgusting name */
ViscoelasticForceTerm *viscoelasticProject_REPViscoelasticAware_veft = NULL;
ViscoelasticRheology *viscoelasticProject_REPViscoelasticAware_rheology = NULL;

const Type Viscoelastic_REPViscoelasticAware_Type = "Viscoelastic_REPViscoelasticAware";

void Viscoelastic_REPViscoelasticAware_AssembleAtParticle(
  RecoveredFeVariable* self,
  ConstitutiveMatrix*  constitutiveMatrix,
  int                  particleIndex,
  IntegrationPoint*    particle,
  int                  lElement_I,
  double*              globalCoord,
  double**             GNx,
  double               detJac,
  double***            Hi_Mat,
  double**             Fi_Mat )
{

  SymmetricTensor p_StrainRate, p_Stress, stressCopy, strainCopy;
  RheologyMaterial* material = (RheologyMaterial*) IntegrationPointsSwarm_GetMaterialOn( (IntegrationPointsSwarm*)constitutiveMatrix->integrationSwarm, particle );
  int dim  = self->dim;
  int dofThatExist = self->fieldComponentCount;
  int order = self->orderOfInterpolation;
  int dof_I;
  double pressure, factor;
  double* pVec   = self->pVec;
  double** pMatrix = self->pMatrix;
  double** tmpC = self->tmpC;

  /* integration factor */
  factor = particle->weight * detJac;

  // calculate polynomial for REP
  self->_makePoly( globalCoord, pVec );

  // get raw strain rate values on particle
  FeVariable_InterpolateWithinElement( self->rawField, lElement_I, particle->xi, p_StrainRate );

  // get the constitutive matrix on the particle
  // and calculate raw deviatoric stress

  assert( constitutiveMatrix );

  if( constitutiveMatrix->storeConstitutiveMatrix )
  {
    // use the stored constMatrix to match current velocity and pressure solution
    ConstitutiveMatrix_GetStoredMatrixOnParticle( constitutiveMatrix, particle, tmpC );
    // perform dev. stress calculation 
    if (dim == 2 ) _TmpConstMat2D_CalculateStress( tmpC, p_StrainRate, p_Stress, constitutiveMatrix->isDiagonal );
    else _TmpConstMat3D_CalculateStress( tmpC, p_StrainRate, p_Stress, constitutiveMatrix->isDiagonal );
  }
  else
  {
    /* build total constitutive matrix - could disagree with current velcoity and pressure solution */
    ConstitutiveMatrix_Assemble( constitutiveMatrix, lElement_I, particleIndex, particle );
    ConstitutiveMatrix_CalculateStress( constitutiveMatrix, p_StrainRate, p_Stress );
  }

  if( material->compressible )
  {
    /* Check if material is compressible, if so adjust direct stress components, here it's assumed the
    compressible material has an isotropic viscosity */
    double du_dx, dv_dy, dw_dz;
    double compressibleBit;
    TensorArray p_velGradField;
    FeVariable_InterpolateWithinElement( self->velGradField, lElement_I, particle->xi, p_velGradField );
    du_dx = p_velGradField[0];

    if( dim == 2 )
    {
      dv_dy = p_velGradField[2];
      compressibleBit = (-2/3)*constitutiveMatrix->matrixData[0][0]*(du_dx+dv_dy);
      p_Stress[0] = p_Stress[0] + compressibleBit;
      p_Stress[1] = p_Stress[1] + compressibleBit;
    }
    else
    {
      dv_dy = p_velGradField[4];
      dw_dz = p_velGradField[8];
      compressibleBit = (-2/3)*constitutiveMatrix->matrixData[0][0]*(du_dx+dv_dy+dw_dz);
      p_Stress[0] = p_Stress[0] + compressibleBit;
      p_Stress[1] = p_Stress[1] + compressibleBit;
      p_Stress[2] = p_Stress[2] + compressibleBit;
    }
  }

  /* THE IMPORTANT BIT */
  {
    /* get viscoelastic particle extension */

    MaterialPointsSwarm* materialPointsSwarm;
    MaterialPoint* materialPoint = OneToOneMapper_GetMaterialPoint( ((IntegrationPointsSwarm*)constitutiveMatrix->integrationSwarm)->mapper, particle, &materialPointsSwarm );

    Viscoelastic_Particle* veExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, viscoelasticProject_REPViscoelasticAware_veft->particleExtHandle );
    /* get viscoelastic rheology */
    ViscoelasticRheology* veRheology = viscoelasticProject_REPViscoelasticAware_rheology;
    double mu = veRheology->mu;
    double dt_e = veRheology->elasticTimeStep;
    double factor = 0;
    if( dim == 2 )
    {
      if (constitutiveMatrix->storeConstitutiveMatrix) factor = tmpC[2][2]/(mu*dt_e);
      else  factor = constitutiveMatrix->matrixData[2][2]/(mu*dt_e);
    }
    else
    {
      if (constitutiveMatrix->storeConstitutiveMatrix) factor = tmpC[5][5]/(mu*dt_e);
      else factor = constitutiveMatrix->matrixData[5][5]/(mu*dt_e);
    }

    if( !self->recoverStrain ) {
      p_Stress[0] = p_Stress[0] + (factor) * veExt->prevStress[0];
      p_Stress[1] = p_Stress[1] + (factor) * veExt->prevStress[1];
      p_Stress[2] = p_Stress[2] + (factor) * veExt->prevStress[2];
      if( dim == 3 )
      {
	p_Stress[3] = p_Stress[3] + (factor) * veExt->prevStress[3];
	p_Stress[4] = p_Stress[4] + (factor) * veExt->prevStress[4];
	p_Stress[5] = p_Stress[5] + (factor) * veExt->prevStress[5];
      }
    }
  }

  /* If pressure exists, then recovered total pressure by taking the raw pressure field as input */
  if ( self->rawPressureField != NULL && !self->recoverStrain )
  {
    FeVariable_InterpolateWithinElement( self->rawPressureField, lElement_I, particle->xi, &pressure );
    p_Stress[0] = p_Stress[0] - pressure; // xx
    p_Stress[1] = p_Stress[1] - pressure; // yy
    if( dim == 3 )
      p_Stress[2] = p_Stress[2] - pressure; // zz
  }

  // make copies of stress and strain tensors
  if( self->recoverStrain ) {
    memcpy( strainCopy, p_StrainRate, dofThatExist*sizeof(double) );
    if( !constitutiveMatrix->storeConstitutiveMatrix )
    {
      /* copy constitutive matrix */
      Index ii;
      for(ii = 0 ; ii < dofThatExist ; ii++)
	memcpy( tmpC[ii], constitutiveMatrix->matrixData[ii], dofThatExist*sizeof(double) );

    }
    if( dim == 2 )
    {
      /* modify stored constitutive matrix to account for symmetry (see CalculateStress routines)  */
      tmpC[0][2] *= 2.0;
      tmpC[1][2] *= 2.0;
      tmpC[2][2] *= 2.0;
    }
    else
    {
      tmpC[0][3] *= 2.0;
      tmpC[1][3] *= 2.0;
      tmpC[2][3] *= 2.0;
      tmpC[3][3] *= 2.0;
      tmpC[4][3] *= 2.0;
      tmpC[5][3] *= 2.0;

      tmpC[0][4] *= 2.0;
      tmpC[1][4] *= 2.0;
      tmpC[2][4] *= 2.0;
      tmpC[3][4] *= 2.0;
      tmpC[4][4] *= 2.0;
      tmpC[5][4] *= 2.0;

      tmpC[0][5] *= 2.0;
      tmpC[1][5] *= 2.0;
      tmpC[2][5] *= 2.0;
      tmpC[3][5] *= 2.0;
      tmpC[4][5] *= 2.0;
      tmpC[5][5] *= 2.0;
    }

  }

  memcpy( stressCopy, p_Stress, dofThatExist*sizeof(double) );

  for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ )
  {
    /* now assemble things per dof, see eq. 20 in
     * B.Boroomand & O.C.Zienkiewicz,
     * "An Improved REP Recovery and the Effectivity Robustness Test",
     * Int. J. for Numerical Methods in Engineering, vol. 40, pages 3247-3277, 1997. */

    /* Setup raw stress tensor for RHS */
    memset( p_Stress, 0, dofThatExist*sizeof(double) );
    p_Stress[dof_I] = stressCopy[dof_I];

    /* Setup polynomial matrix for creating LHS */
    ZeroMatrix( pMatrix, dofThatExist, order );
    memcpy( pMatrix[dof_I], pVec, order*sizeof(double) );

    if( self->recoverStrain )
    {
      /* Multiply in tmpC so that we solve for strainrates */
      NonSquareMatrix_MultiplicationByNonSquareMatrix(
	tmpC, dofThatExist, dofThatExist,
	pMatrix, dofThatExist, order,
	self->CPmat);

      self->_calcHi( self, GNx, self->CPmat, factor, Hi_Mat[dof_I] );
    }
    else
    {
      /* cal Hi, build lefthand Operator */
      self->_calcHi( self, GNx, pMatrix, factor, Hi_Mat[dof_I] );
    }
    /* cal Fi, build righthand Operator */
    self->_calcFi( self, GNx, p_Stress, factor, Fi_Mat[dof_I] );
  }
}


void _Viscoelastic_REPViscoelasticAware_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )
{
  UnderworldContext* context = NULL;
  REP_Algorithm* rep_algorithm = NULL;
  int field_i;

  context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  );
  rep_algorithm = (REP_Algorithm*)Stg_ComponentFactory_ConstructByName( cf, (Name)"REP", REP_Algorithm, True, data  );
  viscoelasticProject_REPViscoelasticAware_veft = Stg_ComponentFactory_ConstructByName( cf, (Name)"viscoelasticForceTerm", ViscoelasticForceTerm, True, data  );
  viscoelasticProject_REPViscoelasticAware_rheology = Stg_ComponentFactory_ConstructByName( cf, (Name)"viscoelasticRheology", ViscoelasticRheology, True, data  );

  /* for each recovered field set it's function ptr to use the
  * Viscoelastic_REPViscoelasticAware_AssembleAtParticle()
  */
  for( field_i = 0 ; field_i < rep_algorithm->repFieldCount ; field_i++ )
  {
    rep_algorithm->repFieldList[field_i]->_assembleOnParticle = Viscoelastic_REPViscoelasticAware_AssembleAtParticle;
  }
}


void* _Viscoelastic_REPViscoelasticAware_DefaultNew( Name name )
{
  return Codelet_New(
           Viscoelastic_REPViscoelasticAware_Type,
           _Viscoelastic_REPViscoelasticAware_DefaultNew,
           _Viscoelastic_REPViscoelasticAware_AssignFromXML,
           _Codelet_Build,
           _Codelet_Initialise,
           _Codelet_Execute,
           _Codelet_Destroy,
           name );
}


Index Viscoelastic_REPViscoelasticAware_Register( PluginsManager* pluginsManager )
{
  return PluginsManager_Submit( pluginsManager, Viscoelastic_REPViscoelasticAware_Type, (Name)"0", _Viscoelastic_REPViscoelasticAware_DefaultNew  );
}


