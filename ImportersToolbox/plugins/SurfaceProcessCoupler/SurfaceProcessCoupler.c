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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <ImportersToolbox/Base/Base.h>

#include <unistd.h>

const Type ImportersToolbox_SurfaceProcessCoupler_Type = "ImportersToolbox_SurfaceProcessCoupler";

typedef struct { 
   __Codelet
   char *sync_folder;            // path to sync folder
   char *maestro_path;           // path to maestro file
   char *ascii_path;             // path to uw ascii file output
   unsigned int sleep_interval;  // sleep interval for wait for SP code
   MaterialPointsSwarm* pts;     // passive tracer swarm
   MaterialPointsSwarm* ms;      // MaterialPointsSwarm
   int vertAxis;                 // vertical axis definition
   double SP_tracer_height;      // height to place SP tracers
   BelowHeightField* belowSurfaceShape;    // representation of the upper surface
   RheologyMaterial* air_material;
   RheologyMaterial* sediment_material;
   SwarmVariable* SP_globalid_var; // the globalid swarm variable
   // remember time information
   double start_time;
   double sync_time;
   double absolute_sync_time;
   int poll_maestro;             //state flag for poll the maestro file
} ImportersToolbox_SurfaceProcessCoupler;

   ImportersToolbox_SurfaceProcessCoupler* global_self=NULL;


void _ImportersToolbox_SurfaceProcessCoupler_load_surface_convert_materials( ImportersToolbox_SurfaceProcessCoupler* self ) {

   /**
     This functions uses the surface process code's surface to modify the materials points.
     
     Material points ABOVE the surface are converted to 'air'
     Material point that were 'air' and are BELOW the surface are converted to 'sediment'.
     */

   MaterialPointsSwarm* ms=NULL;
   MaterialPoint *mp=NULL;
   int ii;
   int local_num_air, local_num_convert;
   int global_num_air, global_num_convert;

   ms = self->ms;
   local_num_air=local_num_convert=0;
   global_num_air=global_num_convert=0;

   assert(ms);

   for( ii=0 ; ii<ms->particleLocalCount; ii++ ) {
      mp = (MaterialPoint*)Swarm_ParticleAt( ms, ii);

      /* is particle outside belowSurfaceShape */ 
      if( !Stg_Shape_IsCoordInside( self->belowSurfaceShape, mp->coord ) ) {
         /* make it air */
         mp->materialIndex = self->air_material->index;
         local_num_air++;
         continue;
      } 

      /* if particle is inside shape, and if also air */ 
      if( mp->materialIndex == self->air_material->index ) {
         /* make it sediment */
         mp->materialIndex = self->sediment_material->index;
         local_num_convert++;
      }
   }

   /* communicate function behaviour */
   (void)MPI_Reduce( &local_num_air, &global_num_air, 1, MPI_INT,
         MPI_SUM, 0, self->context->communicator);
   (void)MPI_Reduce( &local_num_convert, &global_num_convert, 1, MPI_INT,
         MPI_SUM, 0, self->context->communicator);

   if( self->context->rank == 0 ) {
      Journal_Printf( global_info_stream,
            "\nUnderworld converted %d particles to sediment, air population is %d\n",
            global_num_convert, global_num_air );
   }
}

void _ImportersToolbox_SurfaceProcessCoupler_reset_SP_particles( ImportersToolbox_SurfaceProcessCoupler* self ) {
   /**
     Re-initialise the passive tracers representing the surface process's surface underneath uw's top surface
     assumes the swarm only advects vertically for now
     */
   MaterialPointsSwarm* ms=NULL;
   MaterialPoint *mp=NULL;
   int ii;

   ms = self->pts;
   assert(ms);

   /* change vertical height of each local particle */
   for( ii=0 ; ii<ms->particleLocalCount; ii++ ) {
      mp = (MaterialPoint*)Swarm_ParticleAt( ms, ii);
      mp->coord[self->vertAxis] = self->SP_tracer_height;
   }

   /* update particle owners */
   Swarm_UpdateAllParticleOwners( ms );

   if( self->context->rank == 0 ) {
      Journal_Printf( global_info_stream,
            "\nUnderworld has reset the SP code's tracers at height %g\n",
            self->SP_tracer_height );
   }

}

void _ImportersToolbox_SurfaceProcessCoupler_wait_for_SP( void* _self ) {
   ImportersToolbox_SurfaceProcessCoupler* self = (ImportersToolbox_SurfaceProcessCoupler*) _self;

   /* if uw isn't waiting for the SP simulation end function */
   if( !self->poll_maestro )
      return;
   
   if( self->context->rank == 0 ) { 
      FILE *fPtr=NULL;
      fPtr = fopen( self->maestro_path, "r" ); 

      Journal_Printf( global_info_stream, "\nWaiting for Surface Process simulation...\n\n");

      // check if file exists
      Journal_Firewall( fPtr!=NULL, global_debug_stream,
            "Error: in %s. Can't open the maestro file. Path expected is '%s'. Possibly solution is to change the 'sync_folder' parameter",
            __func__, self->sync_folder);
      fclose( fPtr );

      while( ( fPtr = fopen(self->maestro_path, "r") ) ) {
         int c = fgetc( fPtr );
         if( c == 'U' ) {
            fclose( fPtr );
            break;
         }

         fclose( fPtr );
         //wait
         sleep(self->sleep_interval);
      }
   }

   // syncing parallel execution
   MPI_Barrier(self->context->communicator);
   Journal_Printf( global_info_stream, "\nUnderworld is go...\n\n");

   self->poll_maestro=0;

   /* FORCE both the re-initialisation of the SP VoxelField */
   _VoxelFieldVariable_Initialise( self->belowSurfaceShape->heightField, NULL );
   /* if there is a need to re-init the tracers
      _VoxelParticleLayout_InitialiseParticles( self->pts->particleLayout, self->pts ),
      but it errors now */

   _ImportersToolbox_SurfaceProcessCoupler_reset_SP_particles( self );

   /* if not the 1st timestep and not a restart time step we convert things */
   if( self->context->timeStep != 0 &&
         !(self->context->timeStepSinceJobRestart==0 && self->context->loadFromCheckPoint ) ) {
      _ImportersToolbox_SurfaceProcessCoupler_load_surface_convert_materials( self );
   }

}

void _ImportersToolbox_SurfaceProcessCoupler_write_L( ImportersToolbox_SurfaceProcessCoupler* self ) {

   MPI_Barrier( self->context->communicator );

   if( self->context->rank == 0 ) {
      Journal_Printf( global_info_stream,
            "\nUnderworld has written file %s for the Surface Process simulation\n", self->ascii_path);

      FILE *fPtr=NULL;
      fPtr = fopen( self->maestro_path, "w" ); 

      // check if file exists
      Journal_Firewall( fPtr!=NULL, global_debug_stream,
            "Error: in %s. Can't open the maestro file. Path expected is '%s'. Possibly solution is to change the 'sync_folder' parameter",
            __func__, self->sync_folder);

      // put the character
      fputc( 'L', fPtr );

      fclose( fPtr );
   }
}

void _ImportersToolbox_SurfaceProcessCoupler_tracer_output( ImportersToolbox_SurfaceProcessCoupler* self ) {
   /**
     Go though the passive tracer swarm and output the global id and change in height 
     */

   MaterialPointsSwarm* ms=NULL;
   MaterialPoint *mp=NULL;
   
   int ierr, rank, nprocs;
   double globalID;
   int ii, lParticleCount, vertAxis;

   FILE* fPtr=NULL;

   MPI_Status status;
   const int FINISHED_WRITING_TAG = 100;
   int canExecute = 0;

   ms = self->pts;
   lParticleCount = ms->particleLocalCount;
   assert(ms);

   vertAxis = self->vertAxis;
   rank = self->context->rank;
   nprocs = self->context->nproc;

   /* check to see if existing uw_ascii output is there. If so delete */
   if( self->context->rank == 0 ) {
      fPtr = fopen( self->ascii_path, "r" );
      if(fPtr) {
         fclose(fPtr);
         Journal_Firewall( remove(self->ascii_path)==0, global_error_stream,
               "Error in %s: Try to delete file at '%s' but i couldn't\n", __func__, self->ascii_path);
      }
   }

   /* wait for go-ahead from process ranked lower than me, to avoid competition writing to file */
   if ( rank != 0  && canExecute == 0 ) {
      ierr=MPI_Recv( &canExecute, 1, MPI_INT, rank - 1, FINISHED_WRITING_TAG, ms->comm, &status );
   }

      /* open in append mode */
      fPtr = fopen( self->ascii_path, "a" );

      for( ii = 0 ; ii < lParticleCount; ii++ ) {
         mp = (MaterialPoint*)Swarm_ParticleAt( ms, ii );

         /* note we pass globalID as a (double*) because the SwarmVariable_ValueAt can only use that */
         SwarmVariable_ValueAt( self->SP_globalid_var, ii, (&globalID) );

         fprintf( fPtr, "%g %.15g\n", globalID, (mp->coord[vertAxis] - self->SP_tracer_height) );
      }

      /* close file */
      fclose(fPtr);

   /* confirms this processor is finshed */
   canExecute = 1;
   /* send go-ahead from process ranked lower than me, to avoid competition writing to file */
   if ( rank != nprocs - 1 ) {
          MPI_Ssend( &canExecute, 1, MPI_INT, rank + 1, FINISHED_WRITING_TAG, ms->comm );
   }
}

void _ImportersToolbox_SurfaceProcessCoupler_write_for_SP( void* _self, void* data ) {
  // UnderworldContext* context = (ImportersToolbox_SurfaceProcessCoupler*)_self;
   ImportersToolbox_SurfaceProcessCoupler* self= global_self;

   if( (self->context->currentTime-self->absolute_sync_time) > 1e-5*self->context->currentTime ) {
      // freak out because _ImportersToolbox_SurfaceProcessCoupler_dt_limiter didn't perform it's job well 
      Journal_Printf( global_info_stream, "Error in %s. Seems the dt limiter for Underworld screwed up and missed the sync time\n",
            __func__);
      abort();
   }
   if( fabs(self->context->currentTime-self->absolute_sync_time) > 1e-5*self->context->currentTime ) { 
      self->poll_maestro=0;
      return; // keep on running underworld
   }

   self->start_time = self->context->currentTime;
   self->absolute_sync_time = self->start_time + self->sync_time;

   // make the SP tracer output
   _ImportersToolbox_SurfaceProcessCoupler_tracer_output( self );

   // write the 'L' to the maestro file
   _ImportersToolbox_SurfaceProcessCoupler_write_L( self );

   self->poll_maestro=1;

}

void _ImportersToolbox_SurfaceProcessCoupler_Destroy( void* _self, void* data ) {
   ImportersToolbox_SurfaceProcessCoupler* self = (ImportersToolbox_SurfaceProcessCoupler*) _self;

   // free allocated memory
   Memory_Free( self->maestro_path );

   _Codelet_Destroy( _self, data );

}

double _ImportersToolbox_SurfaceProcessCoupler_dt_limiter(void* _self, void* data ) {
   ImportersToolbox_SurfaceProcessCoupler* self = (ImportersToolbox_SurfaceProcessCoupler*) _self;
   double max_interval_for_sync;

   /* set a limite on the max dt the coupling will allow, based on the 'sync_time' */
   max_interval_for_sync = self->absolute_sync_time - self->context->currentTime;

   return max_interval_for_sync;
}

void _ImportersToolbox_SurfaceProcessCoupler_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
   ImportersToolbox_SurfaceProcessCoupler* self = (ImportersToolbox_SurfaceProcessCoupler*) component;
   UnderworldContext* uw_context=NULL;

   uw_context              = (UnderworldContext*)  Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"Context", AbstractContext, True, data  ); 
   self->ms                = (MaterialPointsSwarm*)Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  ); 
   self->pts               = (MaterialPointsSwarm*)Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"UWDeformationTrackingSwarm", MaterialPointsSwarm, True, data  ); 
   self->belowSurfaceShape = (BelowHeightField*)   Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"BelowSurfaceShape", BelowHeightField, True, data  ); 
   self->air_material      = (RheologyMaterial*)   Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"air_material", RheologyMaterial, True, data  ); 
   self->sediment_material = (RheologyMaterial*)   Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"sediment_material", RheologyMaterial, True, data  ); 
   self->sync_folder       = Stg_ComponentFactory_PluginGetString(      cf, self, (Dictionary_Entry_Key)"sync_folder", "./" );
   self->sleep_interval    = Stg_ComponentFactory_PluginGetUnsignedInt( cf, self, (Dictionary_Entry_Key)"sleep_interval", 5 );
   self->sync_time         = Stg_ComponentFactory_PluginGetDouble(      cf, self, (Dictionary_Entry_Key)"sync_time", -1 );
   self->vertAxis          = Stg_ComponentFactory_PluginGetInt(         cf, self, (Dictionary_Entry_Key)"VerticalAxis", 1 );
   self->SP_tracer_height  = Stg_ComponentFactory_PluginGetDouble(      cf, self, (Dictionary_Entry_Key)"SP_tracer_vertical_coord", 0 );
   self->SP_globalid_var   = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Dictionary_Entry_Key)"SP_globalid_SwarmVariable", SwarmVariable, True, data );

   if( self->sync_time < 0 ) {
      Journal_Printf( NULL, "Error in function %s\n. The 'sync_time' input parameters is invalid\n" );
      abort();
   }

   global_self = self;
   self->context = (AbstractContext*)uw_context;
   self->start_time = self->context->currentTime;
   self->absolute_sync_time = self->start_time + self->sync_time;
   self->poll_maestro=1;

   /* create hooks on entry points */
   // hook for waiting for SP 
   EP_AppendClassHook( AbstractContext_GetEntryPoint( self->context, AbstractContext_EP_PreSolveClass),
         _ImportersToolbox_SurfaceProcessCoupler_wait_for_SP, self );

   assert( uw_context->calcDtEP );
   EP_AppendClassHook( uw_context->calcDtEP, _ImportersToolbox_SurfaceProcessCoupler_dt_limiter, self );

   // hook for generating I/O for SP 
   EP_AppendClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass), 
         _ImportersToolbox_SurfaceProcessCoupler_write_for_SP, self);


   /* setup path for maestro file */
   self->maestro_path=Memory_Alloc_Array( char, sizeof(char)*( strlen(self->sync_folder)+8 ), Name_Invalid );
   assert( self->maestro_path );
   sprintf( self->maestro_path, "%s/maestro", self->sync_folder );

   /* setup path for ascii file */
   self->ascii_path=Memory_Alloc_Array( char, sizeof(char)*( strlen(self->sync_folder)+16 ), Name_Invalid );
   assert( self->ascii_path );
   sprintf( self->ascii_path, "%s/uw_output.ascii", self->sync_folder );

}

void* _ImportersToolbox_SurfaceProcessCoupler_DefaultNew( Name name ) {
	/* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof( ImportersToolbox_SurfaceProcessCoupler );
   Type                                                      type = ImportersToolbox_SurfaceProcessCoupler_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ImportersToolbox_SurfaceProcessCoupler_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _ImportersToolbox_SurfaceProcessCoupler_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Codelet_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _ImportersToolbox_SurfaceProcessCoupler_Destroy;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return _Codelet_New(  CODELET_PASSARGS  );
}

Index ImportersToolbox_SurfaceProcessCoupler_Register( PluginsManager* pluginsManager ) {
   return PluginsManager_Submit( pluginsManager, ImportersToolbox_SurfaceProcessCoupler_Type, (Name)"0", _ImportersToolbox_SurfaceProcessCoupler_DefaultNew  );
}


