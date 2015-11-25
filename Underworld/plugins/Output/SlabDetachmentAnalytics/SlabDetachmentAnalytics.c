/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (c) 2005-2010, Monash University, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
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
/** \file
**/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <assert.h>


typedef struct {
   __Codelet
   double             neckingZero;
   GeneralSwarm*      neckingTracers;
   GeneralSwarm*      depthTracers;
} Underworld_SlabDetachmentAnalytics;


const Type Underworld_SlabDetachmentAnalytics_Type = "Underworld_SlabDetachmentAnalytics";
static Bool SlabDetachmentAnalyticsBeenHere;
static Stream* SlabDetachmentAnalyticsStream;

void SlabDetachmentAnalytics( Underworld_SlabDetachmentAnalytics* self, UnderworldContext* context )
{
  static unsigned     ppcIntegralCount;
  static int          ppcIntegralIndex[30];
  unsigned     ii;

  if (!SlabDetachmentAnalyticsBeenHere ) {
    Name                 filename;
    Bool                 fileOpened;
    Stream*              errorStream  = Journal_Register( Error_Type, (Name)CURR_MODULE_NAME  );
    SlabDetachmentAnalyticsStream = Journal_Register( Info_Type, (Name)"SlabDetachmentAnalytics" );

    /* Set up stream */
    Stg_asprintf( &filename, "SlabDetachmentAnalytics.dat" );
    /* Open File */
    if ( context->rank == 0 ) {
      if ( context->loadFromCheckPoint == False ) {
        /* Always overwrite the file if starting a new run */
        fileOpened = Stream_RedirectFile_WithPrependedPath( SlabDetachmentAnalyticsStream, context->outputPath, filename );
      } else {
        /* Just append to the file if doing a restart from checkpoint */
        fileOpened = Stream_AppendFile_WithPrependedPath( SlabDetachmentAnalyticsStream, context->outputPath, filename );
      }
      Journal_Firewall( fileOpened, errorStream,
                        "Could not open file %s/%s. Possibly directory %s does not exist or is not writable.\n"
                        "Check 'outputPath' in input file.\n", context->outputPath, filename, context->outputPath );
    }
    Memory_Free( filename );
    Stream_SetAutoFlush( SlabDetachmentAnalyticsStream, True );
    /* get all PpcIntegral components */
    ppcIntegralCount=0;
    for(ii=0; ii< LiveComponentRegister_GetCount( context->CF->LCRegister ); ii++) {
      if( Stg_Class_IsInstance( LiveComponentRegister_At( context->CF->LCRegister, ii ), PpcIntegral_Type) ) {
        ppcIntegralIndex[ppcIntegralCount] = ii;
        ppcIntegralCount++;
      }
    }

    /* Print header to stream */
    if(context->loadFromCheckPoint == False) {
      Journal_Printf( SlabDetachmentAnalyticsStream,
                      "#       Timestep            Time     NeckingMax");
      for(ii=0; ii< ppcIntegralCount; ii++) {
        Stg_Component* comp = LiveComponentRegister_At( context->CF->LCRegister, ppcIntegralIndex[ii] );
        Journal_Printf( SlabDetachmentAnalyticsStream, "%16s", comp->name );
      }
      Journal_Printf( SlabDetachmentAnalyticsStream, "\n" );
    }
    SlabDetachmentAnalyticsBeenHere = True;
  }

  Journal_Printf( SlabDetachmentAnalyticsStream, "    %12.6g    %12.6g",	(double)context->timeStep, context->currentTime );
  // get vd

  GlobalParticle* particle=NULL;
  GeneralSwarm* swarm=self->neckingTracers;
  double neckingZero = self->neckingZero;
  double localMax[3]={0,0,0};
  double *globalList=NULL;
  double neckingMax[3]={0,0,0};
  int p_i;
  int dim=swarm->dim;
  /* get necking max displacement */
  for( p_i=0; p_i < swarm->particleLocalCount; p_i++ ) {
     particle = (GlobalParticle*)Swarm_ParticleAt( swarm, p_i );

     // check the max. horizontal displacement from original position
     if( (particle->coord[0]-neckingZero) > localMax[0] ) {
        memcpy( localMax, particle->coord, sizeof(double)*dim );
     }
  }

  /** collect all local max necking coord to rank 0 **/

  /* create receive buffer on rank 0 */
  if( context->rank==0) {
     globalList = malloc( sizeof(double)*dim*context->nproc );
  }

  /* communicate max coord to rank 0 */
  (void) MPI_Gather( localMax, dim, MPI_DOUBLE, 
              globalList, dim, MPI_DOUBLE, 
              0, MPI_COMM_WORLD );

  /* search for max horizontal displacement */
  if( context->rank==0) {
     for( p_i=0; p_i < context->nproc; p_i++ ) {
        if( globalList[p_i*dim] > neckingMax[0] ) {
           memcpy( neckingMax, &(globalList[p_i*dim]), sizeof(double) * dim );
        }
     }
    Journal_Printf( SlabDetachmentAnalyticsStream, "    (%12.6g %12.6g)", neckingMax[0], neckingMax[1] );
    free( globalList );
  }




  // get slab depth
  // get depth of necking
  for(ii=0; ii< ppcIntegralCount; ii++) {
    Stg_Component* comp = LiveComponentRegister_At( context->CF->LCRegister, ppcIntegralIndex[ii] );
    Journal_Printf( SlabDetachmentAnalyticsStream, "    %12.6g", PpcIntegral_Integrate(comp) );
  }
  Journal_Printf( SlabDetachmentAnalyticsStream, "\n" );

}


void _Underworld_SlabDetachmentAnalytics_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )
{
  Underworld_SlabDetachmentAnalytics* self = (Underworld_SlabDetachmentAnalytics*)component;
  AbstractContext* context;
  SlabDetachmentAnalyticsBeenHere = False;

  context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", AbstractContext, True, data  );
  self->neckingZero = Stg_ComponentFactory_PluginGetDouble( cf, self, (Name)"NeckingZero", 460e3  );
  self->neckingTracers = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"NeckingTracers", GeneralSwarm, True, data  );
  self->depthTracers = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Name)"DepthTracers", GeneralSwarm, True, data  );

 // ContextEP_Append( context, AbstractContext_EP_FrequentOutput, SlabDetachmentAnalytics );
   EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ),	SlabDetachmentAnalytics, self );

}


void* _Underworld_SlabDetachmentAnalytics_DefaultNew( Name name )
{
  /* Variables set in this function */
  SizeT                                              _sizeOfSelf = sizeof( Underworld_SlabDetachmentAnalytics );
  Type                                                      type = Underworld_SlabDetachmentAnalytics_Type;
  Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
  Stg_Class_PrintFunction*                                _print = _Codelet_Print;
  Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
  Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_SlabDetachmentAnalytics_DefaultNew;
  Stg_Component_ConstructFunction*                    _construct = _Underworld_SlabDetachmentAnalytics_AssignFromXML;
  Stg_Component_BuildFunction*                            _build = _Codelet_Build;
  Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
  Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
  Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

  /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
  AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

  return _Codelet_New(  CODELET_PASSARGS  );
}


Index Underworld_SlabDetachmentAnalytics_Register( PluginsManager* pluginsManager )
{
  Index result;

  result = PluginsManager_Submit( pluginsManager, Underworld_SlabDetachmentAnalytics_Type, (Name)"0", _Underworld_SlabDetachmentAnalytics_DefaultNew  );

  return result;
}


