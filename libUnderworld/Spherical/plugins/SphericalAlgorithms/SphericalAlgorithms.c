/**
	This plugin was created to easily switch out cartesian based algorithms for spherical based algorithms

	The advantages of using a plugin whilst developing the spherical code is that i don't have to spend too much time considering object design as I have access to all data structures in this plugin.
	**/


#include <mpi.h>
#include <float.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Spherical/Spherical.h>

typedef struct
{
   __Codelet
   Bool asbr; // allow solid body rotation, i.e. keep the null space in.
   Mesh *mesh;
   FeVariable *velvar;
   TimeIntegrator *timeIntegrator;
   Stokes_SLE *sle;
   /* THINGS TO SAVE */
} SphericalAlgorithms;

const Type SphericalAlgorithms_Type = "Spherical_SphericalAlgorithms";

/* new algorithms for spherical underworld */
double Mesh_HexType_GetMinimumSeparationGeneral( void* hexType, unsigned elInd, double* perDim );
void _SphericalPeriodicBoundariesManager_Build( void* periodicBCsManager, void* data );
void _SphericalPeriodicBoundariesManager_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I );
void _TimeIntegrator_ExecuteRK2Spherical( void* timeIntegrator, void* data );
void TimeIntegrand_SecondOrderSpherical( void* timeIntegrand, Variable* startValue, double dt );
void _Spherical_RTP_StokesSLE_Initialise( void* component, void* data );
void _SphericalSystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* _context );

void _SphericalAlgorithms_CalibrateVelocity( void* _self ) {
  // SystemLinearEquations_UpdateSolutionOntoNodes();
   Spherical_FeVariable_NonAABCsCalibration( (FeVariable*)_self );
}

void SphericalAlgorithms_SetAlgorithms( SphericalAlgorithms* self )
{
   SphericalGenerator* generator = NULL;
   FiniteElementContext* context=NULL;
   Swarm_Register* sr = Swarm_Register_GetSwarm_Register();
   unsigned int s_i, nSwarms=0;

   self->sle = (Stokes_SLE*)LiveComponentRegister_Get( LiveComponentRegister_GetLiveComponentRegister(), "stokesEqn" );
   self->sle->_updateSolutionOntoNodes = _SphericalSystemLinearEquations_UpdateSolutionOntoNodes;

   context = (FiniteElementContext*)self->velvar->context;

   /* Set the velocity feVariable's calibration function for re-rotating boundary dofs post solve */
   self->velvar->_calibrateBCValues = Spherical_FeVariable_NonAABCsCalibration;

   // is a spherical coordinate mesh check periodic conditions for appropriate swarm advection, ie periodic or not
   if( Stg_Class_IsInstance( self->mesh->generator, SphericalGenerator_Type ) ||
         Stg_Class_IsInstance( self->mesh->generator, Q2SphericalGenerator_Type ) )
   {
      Swarm* swarm=NULL;
      generator = (SphericalGenerator*)self->mesh->generator;

      assert( sr );

      // if there is a material swarm and a periodic mesh change the periodic advection algorithm
      if( generator->periodic[0] || generator->periodic[1] || generator->periodic[2] )
      {
         // get number of swarms
         nSwarms = Swarm_Register_GetCount( sr );

         // for each swarm swap periodic algorithms
         for( s_i = 0 ; s_i < nSwarms; s_i++)
         {
            swarm = Swarm_Register_At( sr, s_i );
            if( Stg_Class_IsInstance( swarm, MaterialPointsSwarm_Type ) && swarm->isAdvecting )
            {
               // if there is no periodic boundaries manager creat it now
               if( ((MaterialPointsSwarm*)swarm)->swarmAdvector->periodicBCsManager == NULL ) {
                  ((MaterialPointsSwarm*)swarm)->swarmAdvector->periodicBCsManager = PeriodicBoundariesManager_New( 
                                             "periodicBCsManager", 
                                             (PICelleratorContext*)context, 
                                             (Mesh*)self->mesh, 
                                             (Swarm*)swarm, NULL );
               }
               ((MaterialPointsSwarm*)swarm)->swarmAdvector->periodicBCsManager->_build = _SphericalPeriodicBoundariesManager_Build;
               ((MaterialPointsSwarm*)swarm)->swarmAdvector->periodicBCsManager->_updateParticle = _SphericalPeriodicBoundariesManager_UpdateParticle;
            }
         }
      }
      Spherical_Get_RotationMatrixIJK = &Spherical_GetRotationMatrixIJK_SphericalNodes;

      // add removal of null space unless we allow it
      if( !self->asbr ) {
         // hackish way of getting the Stokes_SLE. We need it to put the NULL space vector on it
         self->sle->_initialise = _Spherical_RTP_StokesSLE_Initialise;
      }
   }
   if( Stg_Class_IsInstance( self->mesh->generator, ProjectionGenerator_Type ) ||
         Stg_Class_IsInstance( self->mesh->generator, Q2ProjectionGenerator_Type ) || 
         Stg_Class_IsInstance( self->mesh->generator, RSGenerator_Type ) ) 
   {

      Swarm* swarm=NULL;
      generator = (SphericalGenerator*)self->mesh->generator;

      assert( sr );
      // get number of swarms
      nSwarms = Swarm_Register_GetCount( sr );

      // for each swarm swap periodic algorithms
      for( s_i = 0 ; s_i < nSwarms; s_i++)
      {
         swarm = Swarm_Register_At( sr, s_i );
         if( Stg_Class_IsInstance( swarm, MaterialPointsSwarm_Type ) && swarm->isAdvecting )
         {
            ((MaterialPointsSwarm*)swarm)->swarmAdvector->_calculateTimeDeriv = _SwarmAdvector_TimeDeriv_Quicker4IrregularMesh;
         }
      }
      Spherical_Get_RotationMatrixIJK = &Spherical_GetRotationMatrixIJK_ProjectionNodes;
   }

   if( context->maxTimeSteps == -1 ) {
      EntryPoint_PrependClassHook_AlwaysFirst( AbstractContext_GetEntryPoint( context, AbstractContext_EP_DumpClass ), "CalibrateVelocity", _SphericalAlgorithms_CalibrateVelocity, context->type, self->velvar );
   }
}

void _SphericalSystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* _context ) {
   SystemLinearEquations*	self = (SystemLinearEquations*)sle;
   SolutionVector_Index	solnVec_I;
   SolutionVector*		currentSolnVec;

   // call original - which updates the solution from the vector to the dof
   _SystemLinearEquations_UpdateSolutionOntoNodes( sle, _context );

   for ( solnVec_I=0; solnVec_I < self->solutionVectors->count; solnVec_I++ ) {
      currentSolnVec = (SolutionVector*)self->solutionVectors->data[solnVec_I];

      // check if we need to calibrate non axis aligned BC for this solutions vector's feVariable
      if( currentSolnVec->feVariable->nonAABCs ) {
         //currentSolnVec->feVariable->_calibrateBCValues( currentSolnVec->feVariable );
         //SolutionVector_LoadCurrentFeVariableValuesOntoVector( currentSolnVec ); // Annoying i don't like the fact that the FeVariable dof and the petsc vector are duplicates and I have to reload the dofs onto the Vec

         Stg_Class_IsInstance( sle, Stokes_SLE_Type ); // check for this hack

         Vec sol = currentSolnVec->vector;
         Vec n = ((Stokes_SLE*)sle)->null_vector;

         if( n ) {
            // only do this if we have a valid null-space vector on the SLE
            PetscScalar a1, a2, a;

            VecDot(n,sol, &a1);
            VecDot(n,n, &a2);
            a=-a1/a2;
            VecAXPY(sol, a, n);
            VecDot(n,sol, &a1);
            printf("I have subtracted the null space from %s - n*sol is %g\n", currentSolnVec->feVariable->name, a1 );
         }

         // update the feVariable dof memory with petsc vector
         SolutionVector_UpdateSolutionOntoNodes( currentSolnVec );
         // correct the dof's to be in cartesian coordinates only
         currentSolnVec->feVariable->_calibrateBCValues( currentSolnVec->feVariable );
      }
   }

}

void _Spherical_RTP_StokesSLE_Initialise( void* component, void* data ) {

   Stokes_SLE* self = (Stokes_SLE*) component;
   Vec n;
   FeVariable *velField = NULL;
   FeEquationNumber *eq_num = NULL;
   FeMesh *mesh = NULL;
   unsigned localnodes;
   double *xyz, vel[3];
   int ii,eq;

   // call the original initialise function
   _SystemLinearEquations_Initialise(component, data );


   /**************************************
     Following code builds the 2D annulus
     null vector, n = (-pos_y,pos_x)
     ************************************/
   velField = self->uSolnVec->feVariable;
   eq_num = velField->eqNum;
   mesh = velField->feMesh;

   assert( velField );
   assert( eq_num );
   assert( mesh );

   // then build null-space vector
   VecDuplicate( self->uSolnVec->vector, &n);

   localnodes = Mesh_GetLocalSize(mesh, MT_VERTEX);

   //visit each node and build rotation velocity
   for(ii=0; ii<localnodes ; ii++) {
      // get the xyz position
      xyz = Mesh_GetVertex( mesh, ii );

      //build null vector at node
      vel[0] = -1*xyz[1];
      vel[1] = xyz[0];

      // if we are at boudary node assume solid body rotation vec
      if( IndexSet_IsMember( mesh->bndNodeSet, ii ) ) {
         vel[0] = 0;
         vel[1] = 1;
      }
      // set vx in null vector
      eq = eq_num->destinationArray[ii][0];
      VecSetValue(n,eq,vel[0],INSERT_VALUES);

      // set vy in null vector
      eq = eq_num->destinationArray[ii][1];
      VecSetValue(n,eq,vel[1],INSERT_VALUES);
   }

   VecAssemblyBegin( n );
   VecAssemblyEnd( n );

   self->null_vector = n;

}

void _SphericalAlgorithms_Initialise(void* component, void* data)
{
   SphericalAlgorithms* self = (SphericalAlgorithms*)component;

   /* Revisit Spherical advection later */
   if( Stg_Class_IsInstance( self->mesh->generator, SphericalGenerator_Type ) )
   {
      Swarm_Register* sr = Swarm_Register_GetSwarm_Register();
      int nSwarms, s_i;
      Swarm* swarm=NULL;
      if( self->timeIntegrator )
      {
         self->timeIntegrator->_execute = _TimeIntegrator_ExecuteRK2Spherical;
      }

      // get number of swarms
      nSwarms = Swarm_Register_GetCount( sr );
      // for each MaterialPointsSwarm swap search algorithms
      for( s_i = 0 ; s_i < nSwarms; s_i++)
      {
         swarm = Swarm_Register_At( sr, s_i );
         if( Stg_Class_IsInstance( swarm, MaterialPointsSwarm_Type ) )
         {
            // use better search algorithm for linear elements
            ((MaterialPointsSwarm*)swarm)->cellLayout->_cellOf = _ElementCellLayout_CellOf_Irregular;
         }
      }
   }
}

void _SphericalAlgorithms_Build(void* _self, void* data)
{
}


void _SphericalAlgorithms_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )
{
   SphericalAlgorithms* self = (SphericalAlgorithms*)component;

   self->mesh = (Mesh*)Stg_ComponentFactory_ConstructByName( cf, (Name)"elementMesh", Mesh, True, data  );
   self->velvar = (FeVariable*)Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityField", FeVariable, True, data  );
   self->timeIntegrator = (TimeIntegrator*)Stg_ComponentFactory_ConstructByName( cf, (Name)"timeIntegrator", TimeIntegrator, False, data  );
   // allow solid body rotation - asbr
   self->asbr = Stg_ComponentFactory_GetRootDictBool( cf, (Name)"AllowSolidBodyRotation", False  );

   self->sle = NULL;

   /** Create a function at the end of the Assign From XML phase that will swap algorithms **/
   SphericalAlgorithms_SetAlgorithms(self);

}


void* _SphericalAlgorithms_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(SphericalAlgorithms);
   Type                                                      type = SphericalAlgorithms_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SphericalAlgorithms_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _SphericalAlgorithms_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _SphericalAlgorithms_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _SphericalAlgorithms_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

   /* default value NON_GLOBAL */
   AllocationType nameAllocationType = NON_GLOBAL;

   return _Codelet_New( CODELET_PASSARGS );
}


Index Spherical_SphericalAlgorithms_Register( PluginsManager* pluginsManager )
{
   return PluginsManager_Submit( pluginsManager, SphericalAlgorithms_Type, (Name)"0", _SphericalAlgorithms_DefaultNew  );
}


void _SphericalPeriodicBoundariesManager_Build( void* periodicBCsManager, void* data )
{
   PeriodicBoundariesManager* self = (PeriodicBoundariesManager*)periodicBCsManager;
   PeriodicBoundary*          boundary=NULL;

   Stg_Component_Build( self->swarm, data, False );
   Stg_Component_Build( self->mesh, data, False );
   self->size = 4;
   self->boundaries = Memory_Alloc_Array( PeriodicBoundary, self->size, "PeriodicBoundariesManager->boundaries" );

   if ( Stg_Class_IsInstance( self->mesh->generator, SphericalGenerator_Type ) )
   {
      SphericalGenerator* sphericalGenerator = (SphericalGenerator*) self->mesh->generator;
      Dimension_Index dim_I;

      // we only expected dim_I = 0, i.e the boundary walls in theta to be periodic
      for ( dim_I = 0 ; dim_I < self->swarm->dim ; dim_I++ )
      {
         /* Add boundaries straight from mesh generator */
         if ( sphericalGenerator->periodic[ dim_I ] )
         {
            boundary = &self->boundaries[self->count];
            boundary->axis = dim_I;
            // angles are converted from DEGREES to RADIANS and stored in RADIANS
            boundary->minWall = (M_PI/180.0)*sphericalGenerator->crdMin[dim_I];
            boundary->maxWall = (M_PI/180.0)*sphericalGenerator->crdMax[dim_I];
            boundary->particlesUpdatedMinEndCount=0;
            boundary->particlesUpdatedMaxEndCount=0;
            self->count++;
         }
      }
   }
}

void _SphericalPeriodicBoundariesManager_UpdateParticle( void* periodicBCsManager, Particle_Index lParticle_I )
{
   PeriodicBoundariesManager*	self = (PeriodicBoundariesManager*)periodicBCsManager;
   double               difference = 0.0;
   double               rtp[3];
   GlobalParticle*      particle = NULL;
   PeriodicBoundary*    bInfo = NULL;

   /* ASSUMES: periodic boundaries are in the y axis of the code */
   bInfo = &self->boundaries[0];

   particle = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );

   Spherical_XYZ2RTP2D( particle->coord, rtp );

   int newOwningCell = CellLayout_CellOf( self->swarm->cellLayout, particle );

   if( rtp[1] < bInfo->minWall )
   {
      // calculate how much particle is outside by
      difference = bInfo->minWall - rtp[1];

      // project particle to the outside
      rtp[1] = bInfo->maxWall - difference;
      Spherical_RTP2XYZ(rtp, particle->coord );

      //particle->owningCell = CellLayout_CellOf( self->swarm->cellLayout, particle );
      bInfo->particlesUpdatedMinEndCount++; //is it needed?
   }
   else if (rtp[1] > bInfo->maxWall )
   {
      // calculate how much particle is outside by
      difference = rtp[1] - bInfo->maxWall;

      // project particle to the outside
      rtp[1] = bInfo->minWall + difference;
      Spherical_RTP2XYZ(rtp, particle->coord );

      //particle->owningCell = CellLayout_CellOf( self->swarm->cellLayout, particle );
      bInfo->particlesUpdatedMaxEndCount++;
   }

   /* TODO: this is a bit of a hack to print this here using the lParticleI = swarm->total - 1, but its
   the only way I can see given this func is part of the SwarmAdvector intermediate. Should really be a
   function on this class that updates all the particles. -- Main.PatrickSunter 15 May 2006 */
   if ( lParticle_I == (self->swarm->particleLocalCount-1) )
   {
      PeriodicBoundary*	boundary = NULL;
      Index					perB_I;

      Journal_DPrintfL( self->debug, 1, "PeriodicBoundariesManager total particles updated:\n" );
      Stream_Indent( self->debug );

      for ( perB_I = 0; perB_I < self->count; perB_I++ )
      {
         boundary = &self->boundaries[perB_I];

         /*
         Journal_DPrintfL( self->debug, 1, "Periodic Boundary in %c Axis {%.2f,%.2f}: %d min end, %d max end\n",
         	IJKTopology_DimNumToDimLetter[boundary->axis], boundary->minWall, boundary->maxWall,
         	boundary->particlesUpdatedMinEndCount, boundary->particlesUpdatedMaxEndCount );
         	*/
         /* Reset the counters for next time */
         boundary->particlesUpdatedMinEndCount = 0;
         boundary->particlesUpdatedMaxEndCount = 0;
      }
      Stream_UnIndent( self->debug );
   }
}

void TimeIntegrand_SecondOrderSpherical( void* timeIntegrand, Variable* startValue, double dt )
{
   TimeIntegrand*	self           = (TimeIntegrand*)timeIntegrand;
   Variable*       variable       = self->variable;
   double*         arrayDataPtr;
   double*         phi_n, *phi_n_rtp;
   double*         phi_tmp_rtp;
   double*         timeDeriv, *vel_rtp;
   Index           component_I;
   Index           componentCount = *variable->dataTypeCounts;
   Index           array_I;
   Index           arrayCount;
   double          startTime      = TimeIntegrator_GetTime( self->timeIntegrator );
   Bool            successFlag = False;
   Stream*         errorStream = Journal_Register( Error_Type, (Name)self->type  );

   timeDeriv = Memory_Alloc_Array( double, componentCount, "Time Deriv" );
   vel_rtp = Memory_Alloc_Array( double, componentCount, "vel_rtp" );
   phi_n = Memory_Alloc_Array( double, componentCount, "phi_n" );
   phi_n_rtp = Memory_Alloc_Array( double, componentCount, "phi_n" );
   phi_tmp_rtp = Memory_Alloc_Array( double, componentCount, "phi_tmp_rtp" );

   memset( timeDeriv, 0, componentCount * sizeof( double ) );
   memset( vel_rtp, 0, componentCount * sizeof( double ) );
   memset( phi_n, 0, componentCount * sizeof( double ) );
   memset( phi_n_rtp, 0, componentCount * sizeof( double ) );
   memset( phi_tmp_rtp, 0, componentCount * sizeof( double ) );

   /* Update Variables */
   Variable_Update( variable );
   Variable_Update( startValue );
   arrayCount     = variable->arraySize;

   for ( array_I = 0 ; array_I < arrayCount ; array_I++ )
   {
      arrayDataPtr = Variable_GetPtrDouble( variable, array_I );

      TimeIntegrator_SetTime( self->timeIntegrator, startTime );

      /* get phi_n */
      memcpy( phi_n, arrayDataPtr, sizeof( double ) * componentCount );

      /* get v_0 */
      successFlag = _SwarmAdvector_TimeDeriv_Quicker4IrregularMesh( self, array_I, timeDeriv );

      /* convert phi_n and v_0 from cartesian to spherical */
      Spherical_VectorXYZ2RTP( timeDeriv, arrayDataPtr, componentCount, vel_rtp );
      Spherical_XYZ2RTP2D( phi_n, phi_n_rtp );

      Journal_Firewall( True == successFlag, errorStream,
                        "Error - in %s(), for TimeIntegrand \"%s\" of type %s: When trying to find time "
                        "deriv for item %u in step %u, *failed*.\n",
                        __func__, self->name, self->type, array_I, 1 );

      // advect half a timestep in spherical coords
      for ( component_I = 0 ; component_I < componentCount ; component_I++ )
         phi_tmp_rtp[ component_I ] = phi_n_rtp[ component_I ] + 0.5 * dt * vel_rtp[ component_I ];

      // convert location to cartesian
      Spherical_RTP2XYZ( phi_tmp_rtp, arrayDataPtr );

      TimeIntegrand_Intermediate( self, array_I );

      TimeIntegrator_SetTime( self->timeIntegrator, startTime + 0.5 * dt );

      /* get v_n_1/2 in cartesian coords */
      successFlag = _SwarmAdvector_TimeDeriv_Quicker4IrregularMesh( self, array_I, timeDeriv );

      // convert v_n_1/2 from cartesian to spherical
      Spherical_VectorXYZ2RTP( timeDeriv, arrayDataPtr, componentCount, vel_rtp );

      if ( True == successFlag )
      {
         for ( component_I = 0 ; component_I < componentCount ; component_I++ )
         {
            phi_tmp_rtp[ component_I ] = phi_n_rtp[ component_I ] + dt * vel_rtp[ component_I ];
         }

         // convert location to cartesian
         Spherical_RTP2XYZ( phi_tmp_rtp, arrayDataPtr );

         TimeIntegrand_Intermediate( self, array_I );
      }
      else
      {
         Journal_Firewall( True == self->allowFallbackToFirstOrder, errorStream,
                           "Error - in %s(), for TimeIntegrand \"%s\" of type %s: When trying to find time "
                           "deriv for item %u in step %u, *failed*, and self->allowFallbackToFirstOrder "
                           "not enabled.\n", __func__, self->name, self->type, array_I, 2 );

         _TimeIntegrand_RewindToStartAndApplyFirstOrderUpdate( self,
               arrayDataPtr, phi_n, startTime, dt,
               timeDeriv, array_I );
      }
   }

   Memory_Free( timeDeriv );
   Memory_Free(vel_rtp);
   Memory_Free(phi_n);
   Memory_Free(phi_n_rtp);
   Memory_Free(phi_tmp_rtp);
}

void _TimeIntegrator_ExecuteRK2Spherical( void* timeIntegrator, void* data )
{
   TimeIntegrator*	self = (TimeIntegrator*) timeIntegrator;
   AbstractContext*	context = (AbstractContext*) self->context;
   Index					integrand_I;
   Index					integrandCount = TimeIntegrator_GetCount( self );
   double				dt = AbstractContext_Dt( context );
   TimeIntegrand*	integrand;
   double wallTime,tmin,tmax;

   Journal_DPrintf( self->debug, "In %s for %s '%s'\n", __func__, self->type, self->name );


   wallTime = MPI_Wtime();
   TimeIntegrator_Setup( self );

   for ( integrand_I = 0 ; integrand_I < integrandCount ; integrand_I++ )
   {
      integrand = TimeIntegrator_GetByIndex( self, integrand_I );

      TimeIntegrator_SetTime( self, context->currentTime );

      wallTime = MPI_Wtime();
      TimeIntegrand_SecondOrderSpherical( integrand, integrand->variable, dt );

      wallTime = MPI_Wtime()-wallTime;
      MPI_Reduce( &wallTime, &tmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
      MPI_Reduce( &wallTime, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
      Journal_RPrintf(self->info,"\t2nd order: %35s - %9.4f [min] / %9.4f [max] (secs)\n", integrand->name, tmin, tmax);

   }

   TimeIntegrator_Finalise( self );

}
