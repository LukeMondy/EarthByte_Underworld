/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**   Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**   Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**   AuScope - http://www.auscope.org
**   Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**   Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**   Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**   Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**   Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SystemSetup/SystemSetup.h>

#include "types.h"
#include "AdvectionDiffusionSLE.h"
#include "Multicorrector.h"
#include "Residual.h"
#include "Timestep.h"

#include <assert.h>

/* Textual name of this class */
const Type AdvDiffMulticorrector_Type = "AdvDiffMulticorrector";

AdvDiffMulticorrector* AdvDiffMulticorrector_New( 
   Name            name,
   double          gamma,
   Iteration_Index multiCorrectorIterations )
{
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)_AdvDiffMulticorrector_DefaultNew( name );

   AdvDiffMulticorrector_InitAll( self, gamma, multiCorrectorIterations );
   return self;
}

/* Creation implementation / Virtual constructor */
AdvDiffMulticorrector* _AdvDiffMulticorrector_New( ADVDIFFMULTICORRECTOR_DEFARGS ) {
   AdvDiffMulticorrector* self;
   
   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(AdvDiffMulticorrector) );
   self = (AdvDiffMulticorrector*)_SLE_Solver_New( SLE_SOLVER_PASSARGS );
   
   /* Virtual info */
   
   return self;
}

void _AdvDiffMulticorrector_Init( 
   AdvDiffMulticorrector* self, 
   double                 gamma,
   Iteration_Index        multiCorrectorIterations )
{
   self->gamma = gamma;
   self->multiCorrectorIterations = multiCorrectorIterations;
}

void AdvDiffMulticorrector_InitAll( 
   void*           solver,
   double          gamma,
   Iteration_Index multiCorrectorIterations )
{
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*) solver;

   SLE_Solver_InitAll( self, False, 0 );
   _AdvDiffMulticorrector_Init( self, gamma, multiCorrectorIterations );
}

void _AdvDiffMulticorrector_Delete( void* solver ) {
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;

   //FreeObject( self->matrixSolver );
   Stg_KSPDestroy(&self->matrixSolver );

   _SLE_Solver_Delete( self );
}

void _AdvDiffMulticorrector_Print( void* solver, Stream* stream ) {
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;
   
   _SLE_Solver_Print( self, stream );

   Journal_PrintValue( stream, self->gamma );
   Journal_PrintValue( stream, self->multiCorrectorIterations );
}

void* _AdvDiffMulticorrector_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                             _sizeOfSelf = sizeof(AdvDiffMulticorrector);
   Type                                                     type = AdvDiffMulticorrector_Type;
   Stg_Class_DeleteFunction*                             _delete = _AdvDiffMulticorrector_Delete;
   Stg_Class_PrintFunction*                               _print = _AdvDiffMulticorrector_Print;
   Stg_Class_CopyFunction*                                 _copy = NULL;
   Stg_Component_DefaultConstructorFunction* _defaultConstructor = _AdvDiffMulticorrector_DefaultNew;
   Stg_Component_ConstructFunction*                   _construct = _AdvDiffMulticorrector_AssignFromXML;
   Stg_Component_BuildFunction*                           _build = _AdvDiffMulticorrector_Build;
   Stg_Component_InitialiseFunction*                 _initialise = _AdvDiffMulticorrector_Initialise;
   Stg_Component_ExecuteFunction*                       _execute = _AdvDiffMulticorrector_Execute;
   Stg_Component_DestroyFunction*                       _destroy = _AdvDiffMulticorrector_Destroy;
   SLE_Solver_SolverSetupFunction*                  _solverSetup = _AdvDiffMulticorrector_SolverSetup;
   SLE_Solver_SolveFunction*                              _solve = _AdvDiffMulticorrector_Solve;
   SLE_Solver_GetResidualFunc*                      _getResidual = NULL;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*)_AdvDiffMulticorrector_New( ADVDIFFMULTICORRECTOR_PASSARGS );
}

void _AdvDiffMulticorrector_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data ) {
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;
   double                 gamma;
   Iteration_Index        multiCorrectorIterations;

   /* Construct Parent */
   _SLE_Solver_AssignFromXML( self, cf, data );

   gamma = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"gamma", 0.5 );
   multiCorrectorIterations = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"multiCorrectorIterations", 2 );


   /* 'safetyFactor' is the fraction of the advection timestep the solver 
      will try and reach if the system is diffusion dominated
      and multiple diffusion timesteps are being solved per stokesEqn
    */
   self->safetyFactor = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"safetyFactor", 0.75 );

   _AdvDiffMulticorrector_Init( self, gamma, multiCorrectorIterations );

   if( self->matrixSolver == PETSC_NULL ) {
      KSPCreate( MPI_COMM_WORLD, &self->matrixSolver );
   }
}

void _AdvDiffMulticorrector_Build( void* solver, void* data ) {
   AdvDiffMulticorrector* self = Stg_CheckType( solver, AdvDiffMulticorrector );

   _SLE_Solver_Build( self, data );
}

void _AdvDiffMulticorrector_Initialise( void* solver, void* data ) {
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;

   _SLE_Solver_Initialise( self, data );
}

void _AdvDiffMulticorrector_Execute( void* solver, void* data ) {
   AdvDiffMulticorrector* self   = Stg_CheckType( solver, AdvDiffMulticorrector );
   AdvectionDiffusionSLE* sle = Stg_CheckType( data, AdvectionDiffusionSLE );

   /* Test code for doing multiple Energy equation steps if diffusion is quicker
      process than advection. The follow code makes assumption like:
      */
   double local_diff_t = 0.;
   double diff_t = 0.;
   double adv_t = 0.;
   double fraction = 0.;
   double safetyFactor = 1.0;  // this needs to be 1.0, now, given that the timestep is selected before all this (see timestep.c)

   //safetyFactor = self->safetyFactor;
   local_diff_t = AdvectionDiffusionSLE_DiffusiveTimestep( sle );

   (void)MPI_Allreduce( &local_diff_t, &diff_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
   
   adv_t = self->context->dt;  // this will be the advective timestep


   if (adv_t < 1e-7 ) 
     fraction=0;
   else
     fraction = diff_t/adv_t;	   

   if( diff_t>=adv_t || fraction < 1e-5 ) { 
     if( self->context->rank == 0 ) { 
       printf("\nAdvection is dominating - adv: %g\t dif: %g\n", adv_t, diff_t );
     }
     _SLE_Solver_Execute( solver, data );
   } else {
     //adv_t = safetyFactor*adv_t; // for timestep safety
     if( self->context->rank == 0 ) { 
	     printf("\nDiffusion is dominating - running multiple steps to catch advection timestep \n" );
     }
     //double interval_t = diff_t;
     double interval_t = 0.;
     while( interval_t < adv_t ) {
       local_diff_t = AdvectionDiffusionSLE_DiffusiveTimestep( sle );
       (void)MPI_Allreduce( &local_diff_t, &diff_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

       if ( interval_t + diff_t > adv_t ) {
	       sle->currentDt = diff_t + (adv_t - (interval_t + diff_t));
       }
       else {
	       sle->currentDt = diff_t; // first set the deltaT for the solver. Because it needs to know.
       }
       _SLE_Solver_Execute( solver, data );

       interval_t += sle->currentDt;

       if( self->context->rank == 0 ) { 
           printf( "Diffusive Timestep = %g -- Cumulative Timestep %g -- Advection Time - %g\n", sle->currentDt, interval_t, adv_t ); 
       }

       if( fabs(adv_t - interval_t) < 1e-2 && (diff_t / adv_t > 1e-2) ) {
          // Detect if we have hit float errors, where the diff_t was unable to exactly match the adv_t
	  // If the abs diff between the desired timestep and the current is very
	  // small, break out and say close enough. Make sure that it's because 
	  // the diffusive timestep is actually very small compared to what is should be.
	       if( self->context->rank == 0 ) { 
		   printf( "Finished multithermal stepping: fabs(adv - diff) = %g, diff / adv: %g\n", fabs(adv_t - interval_t), diff_t / adv_t ); 
	       }
	  break;
       }
     }
   }
}

void _AdvDiffMulticorrector_Destroy( void* solver, void* data ) {
   _SLE_Solver_Destroy( solver, data );
}

void _AdvDiffMulticorrector_SolverSetup( void* solver, void* data ) {
   // AdvDiffMulticorrector* self   = Stg_CheckType( solver, AdvDiffMulticorrector );
   AdvectionDiffusionSLE* sle = Stg_CheckType( data, AdvectionDiffusionSLE );
   
   __AdvDiffResidualForceTerm_UpdateLocalMemory( sle );

   /* The following is disabled, as it appears the stiffness matrix Mat is destroyed during the 
      SystemLinearEquations_MatrixSetup call below.. this results in no ksp mat being set effectively.
      Instead we call KSPSetOperators just before solve in _AdvDiffMulticorrector_CalculatePhiDot_Implicit */
   /* if ( self->matrixSolver && Stg_Class_IsInstance( sle->massMatrix, StiffnessMatrix_Type ) ) {
      StiffnessMatrix* massMatrix = Stg_CheckType( sle->massMatrix, StiffnessMatrix );
      KSPSetOperators( self->matrixSolver, massMatrix->matrix, massMatrix->matrix, DIFFERENT_NONZERO_PATTERN );
   } */ 
}

/* See Brooks, Hughes 1982 Section 4.2 
 * All equations refer to this paper if not otherwise indicated */
void _AdvDiffMulticorrector_Solve( void* solver, void* _sle ) {
   AdvDiffMulticorrector* self = (AdvDiffMulticorrector*) solver;
   AdvectionDiffusionSLE* sle = (AdvectionDiffusionSLE*)_sle;
   double                 dt = sle->currentDt;
   Index                  iteration_I;
   Vec                    deltaPhiDot;

   Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

   /* First apply BC's */
   FeVariable_ApplyBCs( sle->phiVector->feVariable, self->context );
   FeVariable_ApplyBCs( sle->phiDotVector->feVariable, self->context );

   /* Put mesh data onto vectors */
   SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->phiVector );
   SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->phiDotVector );

   /* Solve for predictor step */
   AdvDiffMulticorrector_Predictors( self, sle, dt );

   /* Allocate Memory For Corrector Step */
   //Vector_Duplicate( sle->phiVector->vector, (void**)&deltaPhiDot );
   //Vector_SetLocalSize( deltaPhiDot, Vector_GetLocalSize( sle->phiVector->vector ) );
   VecDuplicate( sle->phiVector->vector, &deltaPhiDot );

   /* Multi-corrector Steps */
   for ( iteration_I = 0 ; iteration_I < self->multiCorrectorIterations ; iteration_I++ ) {
      AdvDiffMulticorrector_Solution( self, sle, deltaPhiDot );
      AdvDiffMulticorrector_Correctors( self, sle, deltaPhiDot, dt );

      /* Put solutions onto meshes */
      SolutionVector_UpdateSolutionOntoNodes( sle->phiVector );
      SolutionVector_UpdateSolutionOntoNodes( sle->phiDotVector );

      SystemLinearEquations_ZeroAllVectors( sle, NULL );
   }

   /* Clean Up */
   //FreeObject( deltaPhiDot );
   Stg_VecDestroy(&deltaPhiDot );
}

void ViewPETScVector( Vec vec, Stream* stream ) {
   PetscInt     size;
   PetscScalar* array;
   unsigned     entry_i;

   if( !stream )
      stream = Journal_Register( Info_Type, (Name)"tmp" );

   VecGetLocalSize( vec, &size );
   VecGetArray( vec, &array );
   
   for( entry_i = 0; entry_i < size; entry_i++ )
      Journal_Printf( stream, "\t%u: \t %.12g\n", entry_i, array[entry_i] );

   VecRestoreArray( vec, &array );
}

/** See Eqns. 4.2.3-4 */
void AdvDiffMulticorrector_Predictors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, double dt ) {
   double  factor = dt * ( 1.0 - self->gamma );
   Stream* debugStream = sle->debug;

   Journal_DPrintf( debugStream, "In func %s:\n", __func__ );

   #if DEBUG
   if ( Stream_IsPrintableLevel( debugStream, 3 ) ) {
      Journal_DPrintf( debugStream, "At start of %s:\n", __func__ );
      Stream_Indent( debugStream );

      Journal_PrintValue( debugStream, dt );
      Journal_PrintValue( debugStream, self->gamma );
      Journal_PrintValue( debugStream, factor );

      Journal_DPrintf( debugStream, "Phi:\n" );
      ViewPETScVector( sle->phiVector->vector, debugStream );
      Journal_DPrintf( debugStream, "Phi Dot:\n" );
      ViewPETScVector( sle->phiDotVector->vector, debugStream );

      Stream_UnIndent( debugStream );
   }
   #endif

   /* Calculate Predictor for \phi - 
    * Eq. 4.2.3: \phi_{n+1}^{(0)} = \phi_n + \Delta t(1 - \gamma)\dot \phi_n */
   //Vector_AddScaled( sle->phiVector->vector, factor, sle->phiDotVector->vector ); 
   VecAXPY( sle->phiVector->vector, factor, sle->phiDotVector->vector ); 
   
   /* Calculate Predictor for \dot \phi - 
    * Eq. 4.2.4: \dot \phi_{n+1}^{(0)} = 0 */
   //Vector_Zero( sle->phiDotVector->vector );
   VecSet( sle->phiDotVector->vector, 0.0 );

   #if DEBUG
   if ( Stream_IsPrintableLevel( debugStream, 3 ) ) {
      Journal_DPrintf( debugStream, "At end of %s: Phi is:\n", __func__ );
      ViewPETScVector( sle->phiVector->vector, debugStream );
   }
   #endif
}
   
void AdvDiffMulticorrector_Solution( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) {
   Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

   /* Calculate Residual - See Eq. 4.2.6 */
   SystemLinearEquations_VectorSetup( sle, NULL );
   SystemLinearEquations_MatrixSetup( sle, NULL );

   /* Calculate Mass Matrix out of three options - fully explicit, fully implicit, split operators */
   AdvDiffMulticorrector_CalculatePhiDot( self, sle, deltaPhiDot );

   #if DEBUG
   if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
      Journal_DPrintf( self->debug, "Delta Phi Dot is:\n" );
      ViewPETScVector( deltaPhiDot, self->debug );
   }
   #endif
}

/* Correct \phi and \dot \phi - See Eqns. 4.2.7-8 */
void AdvDiffMulticorrector_Correctors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot, double dt ) {
   double factor = dt * self->gamma;
   
   Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

   /* Add correction to \phi - Eq. 4.2.7 */
   //Vector_AddScaled( sle->phiVector->vector, factor, deltaPhiDot );
   VecAXPY( sle->phiVector->vector, factor, deltaPhiDot );
   
   /* Add correction to \dot \phi - Eq. 4.2.8 */
   //Vector_AddScaled( sle->phiDotVector->vector, 1.0, deltaPhiDot );
   VecAXPY( sle->phiDotVector->vector, 1.0, deltaPhiDot );
}

void AdvDiffMulticorrector_CalculatePhiDot( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) {
   Stg_Component* massMatrix = sle->massMatrix;

   if ( Stg_Class_IsInstance( massMatrix, ForceVector_Type ) ) 
      _AdvDiffMulticorrector_CalculatePhiDot_Explicit( self, sle, deltaPhiDot );
   else if ( Stg_Class_IsInstance( massMatrix, StiffnessMatrix_Type ) )
      _AdvDiffMulticorrector_CalculatePhiDot_Implicit( self, sle, deltaPhiDot );
   else {
      Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type ),
         "Error in func '%s': Cannot understand type '%s' for mass matrix '%s'.\n",
         __func__, massMatrix->name, massMatrix->type );
   }
}

/* Lump all things onto diagonal of matrix - which is stored as a vector - Eq. 4.2.11 */
void _AdvDiffMulticorrector_CalculatePhiDot_Explicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) {
   ForceVector* massMatrix = Stg_CheckType( sle->massMatrix, ForceVector );

   /* Calculate change in \dot \phi - See Eq. 4.2.5 */
   VecPointwiseDivide( deltaPhiDot, sle->residual->vector, massMatrix->vector );
}

void _AdvDiffMulticorrector_CalculatePhiDot_Implicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) {
   Stg_KSPSetOperators( self->matrixSolver, ((StiffnessMatrix*)sle->massMatrix)->matrix, ((StiffnessMatrix*)sle->massMatrix)->matrix, DIFFERENT_NONZERO_PATTERN );
   KSPSolve( self->matrixSolver, sle->residual->vector, deltaPhiDot );
}


