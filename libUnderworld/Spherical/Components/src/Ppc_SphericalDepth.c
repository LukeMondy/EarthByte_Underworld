#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Components/Components.h>

#include <string.h>
#include <math.h>

#include "types.h"
#include "SphericalUtils.h"
#include "Ppc_SphericalDepth.h"


/* Textual name of this class */
const Type Ppc_SphericalDepth_Type = "Ppc_SphericalDepth";

void _Ppc_SphericalDepth_ReferenceHeight( Ppc_SphericalDepth* self, void* data)
{
   /* This function finds the particle's height and communicate to the other procs */
   double lbuffer, gbuffer;
   int ierr;
   if ( self->ms->particleLocalCount == 1 )
   {
      GlobalParticle *mp;
      double rtp[3];
      mp = (GlobalParticle*)Swarm_ParticleAt(self->ms, 0 );
      (self->ms->dim == 2) ? 
         Spherical_XYZ2RTP2D( mp->coord, rtp ) :
         Spherical_XYZ2RTP3D( mp->coord, rtp ); 
      lbuffer = rtp[0];

   }
   else
   {
      lbuffer = 0.0;
   }

   /* Communicate the height, note MPI_SUM is used */
   ierr = MPI_Allreduce( &lbuffer, &gbuffer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   self->refHeight = gbuffer;

   assert( isnan(self->refHeight) == 0 );
}

void _Ppc_SphericalDepth_Init( Ppc_SphericalDepth* self, FeMesh *mesh, Swarm *swarm, double refHeight )
{
   self->mesh = mesh;
   self->ms = swarm;

   if ( self->ms )
   {
      /* for passiveTraverSwarm as reference height */
      EP_PrependClassHook( Context_GetEntryPoint( self->context, AbstractContext_EP_UpdateClass ),
                           _Ppc_SphericalDepth_ReferenceHeight,
                           self );
   }
   else
   {
      self->refHeight = refHeight;
   }
}


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Ppc_SphericalDepth* _Ppc_SphericalDepth_New(  PPC_CONSTANT_DEFARGS  )
{
   Ppc_SphericalDepth* self;

   assert( _sizeOfSelf >= sizeof(Ppc_SphericalDepth) );
   nameAllocationType = NON_GLOBAL;
   self = (Ppc_SphericalDepth*) _Ppc_New(  PPC_PASSARGS  );
   self->_get = _get;
   return self;
}


void* _Ppc_SphericalDepth_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Ppc_SphericalDepth);
   Type                                                      type = Ppc_SphericalDepth_Type;
   Stg_Class_DeleteFunction*                              _delete = _Ppc_SphericalDepth_Delete;
   Stg_Class_PrintFunction*                                _print = _Ppc_SphericalDepth_Print;
   Stg_Class_CopyFunction*                                  _copy = NULL;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Ppc_SphericalDepth_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Ppc_SphericalDepth_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Ppc_SphericalDepth_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Ppc_SphericalDepth_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Ppc_SphericalDepth_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Ppc_SphericalDepth_Destroy;
   AllocationType                              nameAllocationType = NON_GLOBAL;
   Ppc_GetFunction*                                          _get = _Ppc_SphericalDepth_Get;

   return (void*)_Ppc_SphericalDepth_New(  PPC_CONSTANT_PASSARGS  );
}


void _Ppc_SphericalDepth_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;
   FeMesh* mesh = NULL;
   Swarm* ms = NULL;
   double refHeight = 0;;
   Dictionary*	theDictionary;

   /* Construct parent */
   _Ppc_AssignFromXML( self, cf, data );

   /* The dictionary */
   theDictionary = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, (Dictionary_Entry_Key)self->name )  );
   mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Mesh", FeMesh, False, data );
   if ( !mesh )
      mesh = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_linearMesh", FeMesh, False, data );
   if ( !mesh )
      mesh = self->manager->mesh;


   ms = (Swarm*)Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ReferenceSwarm", MaterialPointsSwarm, False, data );

   if ( ms == NULL )
   {
      Bool existRefHeight = Stg_ComponentFactory_TryDouble( cf, self->name, (Dictionary_Entry_Key)"ReferenceHeight", &refHeight );
      if ( existRefHeight == False )
      {
         Journal_Firewall(0, self->error_stream,
                          "\n\nError in %s, need <param name=\"ReferenceHeight\"></param> explicitly defined in %s\n\n",
                          __func__, self->name );
      }
   }

   /* Init */
   _Ppc_SphericalDepth_Init(
      self,
      mesh,
      ms,
      refHeight);
}


void _Ppc_SphericalDepth_Build( void* _self, void* data )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Build parent */
   _Ppc_Build( self, data );

   if ( self->ms ) Stg_Component_Build(self->ms, data, False );
}

void _Ppc_SphericalDepth_Initialise( void* _self, void* data )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Initialize parent */
   _Ppc_Initialise( self, data );

   if ( self->ms )
   {
      Stg_Component_Initialise(self->ms, data, False );
      _Ppc_SphericalDepth_ReferenceHeight( self, NULL );
   }
}

void _Ppc_SphericalDepth_Delete( void* _self )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Delete parent */
   _Ppc_Delete( self );
}

void _Ppc_SphericalDepth_Print( void* _self, Stream* stream )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Print parent */
   _Ppc_Print( self, stream );
}

void _Ppc_SphericalDepth_Execute( void* _self, void* data )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Execute parent */
   _Ppc_Execute( self, data );
}

void _Ppc_SphericalDepth_Destroy( void* _self, void* data )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*)_self;

   /* Destroy parent */
   _Ppc_Destroy( self, data );
}

/*
 * Public functions
 *
 */

int _Ppc_SphericalDepth_Get( void* _self, unsigned lElement_I, IntegrationPoint* particle, double* result )
{
   Ppc_SphericalDepth* self = (Ppc_SphericalDepth*) _self;
   MaterialPoint *mp = NULL;
   MaterialPointsSwarm *ms=NULL;
   double rtp[3];

   if( self->ms ) {
      mp = OneToOneMapper_GetMaterialPoint( self->manager->integrationSwarm->mapper, particle, &ms );
      (self->ms->dim == 2) ? 
         Spherical_XYZ2RTP2D( mp->coord, rtp ) :
         Spherical_XYZ2RTP3D( mp->coord, rtp ); 
   } else {
      double xyz[3];
      FeMesh_CoordLocalToGlobal( self->manager->mesh, lElement_I, particle->xi, xyz );
      (self->manager->integrationSwarm->dim == 2) ?  
         Spherical_XYZ2RTP2D( xyz, rtp ) : 
         Spherical_XYZ2RTP3D( xyz, rtp ) ; 
   }

   result[0] =  self->refHeight - rtp[0];

   return 0;
}
