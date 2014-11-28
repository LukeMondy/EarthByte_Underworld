#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Components.h"

#include <petscblaslapack.h>

/* Textual name of this class */
const Type SphericalForceVector_Type = "SphericalForceVector";

/** Name of this class' entry points */
static const char	ForceVector_assembleForceVectorStr[] = "assembleForceVector";

void* _SphericalForceVector_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(SphericalForceVector);
   Type                                                      type = ForceVector_Type;
   Stg_Class_DeleteFunction*                              _delete = _ForceVector_Delete;
   Stg_Class_PrintFunction*                                _print = _ForceVector_Print;
   Stg_Class_CopyFunction*                                  _copy = _ForceVector_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SphericalForceVector_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _SphericalForceVector_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _SphericalForceVector_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _ForceVector_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _ForceVector_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _ForceVector_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;

   return _SphericalForceVector_New(  SPHERICALFORCEVECTOR_PASSARGS  );
}

SphericalForceVector* _SphericalForceVector_New(  SPHERICALFORCEVECTOR_DEFARGS  )
{
   SphericalForceVector* self;

   /* Allocate memory */
   assert( _sizeOfSelf >= sizeof(SphericalForceVector) );
   self = (SphericalForceVector*)_ForceVector_New(  FORCEVECTOR_PASSARGS  );

   return self;
}


void _SphericalForceVector_Init( void* forceVector )
{
   SphericalForceVector* self = (SphericalForceVector*)  forceVector;

   self->rotMat = NULL;
   self->tmpMat = NULL;

   /* Add default hook to assembleForceVector entry point */
   EP_ReplaceAll( self->assembleForceVector, SphericalForceVector_GlobalAssembly_General );
}

void _SphericalForceVector_Delete( void* forceVector )
{
   ForceVector* self = (ForceVector*)forceVector;

   Journal_DPrintf( self->debug, "In %s - for %s\n", __func__, self->name );

   /* Stg_Class_Delete parent*/
   _ForceVector_Delete( self );
}

void _SphericalForceVector_AssignFromXML( void* forceVector, Stg_ComponentFactory* cf, void* data )
{
   ForceVector*    self = (ForceVector*)forceVector;

   _ForceVector_AssignFromXML( self, cf, data );

   _SphericalForceVector_Init( self );
}

void _SphericalForceVector_Build( void* forceVector, void* data )
{
   SphericalForceVector* self = (SphericalForceVector*)forceVector;

   _ForceVector_Build( self, data );

   // alloc some memory to these guys
   self->rotMat = ReallocArray( self->rotMat, double, (24*24) );
   self->tmpMat = ReallocArray( self->tmpMat, double, (24*24) );
}
#if 0
void _SphericalForceVector_Initialise( void* forceVector, void* data )
{
   SphericalForceVector* self = (SphericalForceVector*)forceVector;

   _ForceVector_Initialise( self, data );
}

#endif
void _SphericalForceVector_Destroy( void* forceVector, void* data )
{
   SphericalForceVector* self = (SphericalForceVector*)forceVector;

   FreeArray( self->rotMat );
   FreeArray( self->tmpMat );
   _ForceVector_Destroy( self, data );
}

void SphericalForceVector_GlobalAssembly_General( void* forceVector )
{
   SphericalForceVector*   self                 = (SphericalForceVector*) forceVector;
   FeVariable*             feVar                = self->feVariable;
   Element_LocalIndex      element_lI;
   Element_LocalIndex      elementLocalCount;
   Node_ElementLocalIndex  nodeCountCurrElement = 0;
   int                     *nodeIdsInCurrElement = NULL;
   Dof_Index               totalDofsThisElement = 0;
   Dof_Index               totalDofsPrevElement = 0;
   Dof_Index               dofCountLastNode     = 0;
   Dof_EquationNumber**    elementLM            = NULL;
   double*                 elForceVecToAdd      = NULL;
   /* For output printing */
   double                  outputPercentage=10;	/* Controls how often to give a status update of assembly progress */
   int                     outputInterval;

   Journal_DPrintf( self->debug, "In %s - for vector \"%s\"\n", __func__, self->name );

   Stream_IndentBranch( StgFEM_Debug );

   if ( Stg_ObjectList_Count( self->forceTermList ) > 0 )
   {
      elementLocalCount = FeMesh_GetElementLocalSize( feVar->feMesh );

      /* Initialise Vector */
      outputInterval = (int)( (outputPercentage/100.0)*(double)(elementLocalCount) );
      if( outputInterval == 0 )
      {
         outputInterval = elementLocalCount;
      }

      for( element_lI = 0; element_lI < elementLocalCount; element_lI++ )
      {
         int nInc, *inc;

         FeMesh_GetElementNodes( feVar->feMesh, element_lI, self->inc );
         nInc = IArray_GetSize( self->inc );
         inc = IArray_GetPtr( self->inc );
         nodeCountCurrElement = nInc;
         /* Get the local node ids */
         nodeIdsInCurrElement = inc;

         /* Set value of elementLM: will automatically just index into global LM table if built */
         elementLM = FeEquationNumber_BuildOneElementLocationMatrix( feVar->eqNum, element_lI );

         /* work out number of dofs at the node, using LM */
         /* Since: Number of entries in LM table for this element = (by defn.) Number of dofs this element */
         dofCountLastNode = feVar->dofLayout->dofCounts[nodeIdsInCurrElement[nodeCountCurrElement-1]];
         totalDofsThisElement = self->totalDofsThisElement = &elementLM[nodeCountCurrElement-1][dofCountLastNode-1] - &elementLM[0][0] + 1;

         if ( totalDofsThisElement > totalDofsPrevElement )
         {
            if (elForceVecToAdd) Memory_Free( elForceVecToAdd );
            Journal_DPrintfL( self->debug, 2, "Reallocating elForceVecToAdd to size %d\n", totalDofsThisElement );
            elForceVecToAdd = Memory_Alloc_Array( double, totalDofsThisElement, "elForceVecToAdd" );
         }

         /* Initialise Values to Zero */
         memset( elForceVecToAdd, 0, totalDofsThisElement * sizeof(double) );

         /* Assemble this element's element force vector: going through each force term in list */
         SphericalForceVector_AssembleElement( self, element_lI, elForceVecToAdd );

         /* When keeping BCs in we come across a bit of a problem in parallel. We're not
            allowed to add entries to the force vector here and then clobber it later with
            an insert in order to set the BC. So, what we'll do is just add zero here, that
            way later we can add the BC and it will be the same as inserting it.
            --- Luke, 20 May 2008 */
         if( !self->feVariable->eqNum->removeBCs )
         {
            DofLayout* dofs;
            int nDofs, curInd;
            int ii, jj;

            dofs = self->feVariable->dofLayout; /* shortcut to the dof layout */
            curInd = 0; /* need a counter to track where we are in the element force vector */
            for( ii = 0; ii < nodeCountCurrElement; ii++ )
            {
               nDofs = dofs->dofCounts[inc[ii]]; /* number of dofs on this node */
               for( jj = 0; jj < nDofs; jj++ )
               {
                  if( !FeVariable_IsBC( self->feVariable, inc[ii], jj ) )
                  {
                     curInd++;
                     continue; /* only need to clear it if it's a bc */
                  }
                  elForceVecToAdd[curInd] = 0.0;
                  curInd++;
               }
            }
         }

         /* Ok, assemble into global matrix */
         VecSetValues( self->vector, totalDofsThisElement, (PetscInt*)elementLM[0], elForceVecToAdd, ADD_VALUES );

#if DEBUG
         if( element_lI % outputInterval == 0 )
         {
            Journal_DPrintfL( self->debug, 2, "done %d percent of global force vector assembly (general) \n",
                              (int)(100.0*((double)element_lI/(double)elementLocalCount)) );
         }
#endif

         /* Cleanup: If we haven't built the big LM for all elements, free the temporary one */
         if ( False == feVar->eqNum->locationMatrixBuilt )
         {
            Memory_Free( elementLM );
         }
         totalDofsPrevElement = totalDofsThisElement;
      }

      Memory_Free( elForceVecToAdd );
   }
   else
   {
      Journal_DPrintf( self->debug, "No ForceTerms registered - returning.\n" );
   }

   /* If we're keeping BCs, insert them into the force vector. */
   if( !feVar->eqNum->removeBCs )
      Assembler_LoopVector( self->bcAsm );

   VecAssemblyBegin( self->vector );
   VecAssemblyEnd( self->vector );

   Stream_UnIndentBranch( StgFEM_Debug );
}

void SphericalForceVector_AssembleElement( void* forceVector, Element_LocalIndex element_lI, double* elForceVecToAdd )
{
   SphericalForceVector   *self                 = (SphericalForceVector*) forceVector;
   FeVariable             *feVar                = self->feVariable;
   Index                   forceTermCount       = Stg_ObjectList_Count( self->forceTermList );
   Index                   forceTerm_I;
   ForceTerm*              forceTerm;
   int ii;

   for ( forceTerm_I = 0 ; forceTerm_I < forceTermCount ; forceTerm_I++ )
   {
      forceTerm = (ForceTerm*) Stg_ObjectList_At( self->forceTermList, forceTerm_I );

      ForceTerm_AssembleElement( forceTerm, (ForceVector*)self, element_lI, elForceVecToAdd );
   }

   if( feVar->nonAABCs && IndexSet_IsMember( feVar->feMesh->bndElementSet, element_lI ) )
   {
      /*
         Perform [tmp] = [R]^T * [elVecToAdd]
         but do it with BLAS (fortran column major ordered) memory layout
         therefore compute: [tmp]^T = [elVecToAdd]^T * [[R]^T]^T
      */

      double* R = self->rotMat;
      double* tmp = self->tmpMat;

      int rowA = self->totalDofsThisElement; // rows in [R]
      int colA = rowA;              // cols in [R]
      int colB = 1;    // cols in [elStiffMatToAdd]

      PetscScalar one=1.0;
      PetscScalar zero=0.0;
      char t='T';
      char n='N';

      // initialise [R] and [tmp] memory
      memset(R,0,rowA*rowA*sizeof(double));
      memset(tmp,0,rowA*rowA*sizeof(double));

      // evaluate [R] for element_lI
      ( self->dim == 3 ) ?
      SphericalStiffnessMatrix_EvaluateRotationMatrix3D( feVar, element_lI, R ):
      SphericalStiffnessMatrix_EvaluateRotationMatrix2D( feVar, element_lI, R );

      //blasMatrixMult( RT, elForceVecToAdd, 24, 1, 24, tmp );
      // [tmp]^T = [elStiffMatToAdd]^T * [R]
      BLASgemm_( &n, &t,
            &colB, &rowA, &colA, 
            &one, elForceVecToAdd, &colB, 
            R, &rowA, 
            &zero, tmp, &colB );
      // copy result into returned memory segment
      memcpy( elForceVecToAdd, tmp, rowA*colB*sizeof(double) );
     
   }
}
