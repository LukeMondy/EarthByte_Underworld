/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005-2010, Monash University 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**       * Redistributions of source code must retain the above copyright notice, 
**          this list of conditions and the following disclaimer.
**       * Redistributions in binary form must reproduce the above copyright 
**         notice, this list of conditions and the following disclaimer in the 
**         documentation and/or other materials provided with the distribution.
**       * Neither the name of the Monash University nor the names of its contributors 
**         may be used to endorse or promote products derived from this software 
**         without specific prior written permission.
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
/** \file 
\details Calculates the velocity root mean square.
**/
#include <mpi.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "SphericalVrms.h"

const Type Underworld_SphericalVrms_Type = "Spherical_SphericalVrms";

void _Underworld_SphericalVrms_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
   Underworld_SphericalVrms* self = (Underworld_SphericalVrms*)component;

   self->context = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"Context", UnderworldContext, True, data );

   self->gaussSwarm = Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"GaussSwarm", Swarm, True, data );

   self->velocityField = Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"VelocityField", FeVariable, True, data );

   Underworld_SphericalVrms_PrintHeaderToFile( self->context );
   ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Underworld_SphericalVrms_Dump );
}

void _Underworld_SphericalVrms_Build( void* component, void* data ) {
   Underworld_SphericalVrms* self = (Underworld_SphericalVrms*)component;

   assert( self );

   Stg_Component_Build( self->gaussSwarm, data, False );
   Stg_Component_Build( self->velocityField, data, False );
   
   _Codelet_Build( self, data );
}

void _Underworld_SphericalVrms_Initialise( void* component, void* data ) {
   Underworld_SphericalVrms* self = (Underworld_SphericalVrms*)component;

   assert( self );

   Stg_Component_Initialise( self->gaussSwarm, data, False );
   Stg_Component_Initialise( self->velocityField, data, False );
   
   _Codelet_Initialise( self, data );
}

void _Underworld_SphericalVrms_Destroy( void* component, void* data ) {
   Underworld_SphericalVrms* self = (Underworld_SphericalVrms*)component;

   assert( self );

   _Codelet_Destroy( self, data );
   
   Stg_Component_Destroy( self->gaussSwarm, data, False );
   Stg_Component_Destroy( self->velocityField, data, False );
}

void* _Underworld_SphericalVrms_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Underworld_SphericalVrms);
   Type                                                      type = Underworld_SphericalVrms_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_SphericalVrms_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Underworld_SphericalVrms_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Underworld_SphericalVrms_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Underworld_SphericalVrms_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Underworld_SphericalVrms_Destroy;

   /* 
    * Variables that are set to ZERO are variables that will be set either by the
    * current _New function or another parent _New function further up the hierachy.
    */

    /* default value NON_GLOBAL */
   AllocationType nameAllocationType = NON_GLOBAL;

   return _Codelet_New( CODELET_PASSARGS );
}

Index Spherical_SphericalVrms_Register( PluginsManager* pluginsManager ) {
   return PluginsManager_Submit( pluginsManager, Underworld_SphericalVrms_Type, (Name)"0", _Underworld_SphericalVrms_DefaultNew );
}

/* Integrate Every Step and dump to file */
void Underworld_SphericalVrms_Dump( void* _context ) {
   UnderworldContext* context = (UnderworldContext* ) _context;
   FeMesh*              mesh;
   FeVariable* velocityField=NULL;
   double             vec[3], rtp[3], *xyz, vrtp[3];
   double             area, gArea, magVrms;
   double             vrms[2], gVrms[2]; // 1st component is radial_vrms, 2nd angular_vrms 
   Dimension_Index    dim = context->dim;
   int ii, nLocalNodes;

   Underworld_SphericalVrms* self;

   self = (Underworld_SphericalVrms*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Underworld_SphericalVrms_Type );
   assert( self );
   velocityField = self->velocityField; assert( velocityField );
   mesh = velocityField->feMesh; assert( mesh );
   nLocalNodes = FeMesh_GetNodeLocalSize( mesh );
   
   // zero integrals
   vrms[0]=vrms[1]=area=0;

   for( ii=0; ii<nLocalNodes; ii++ ) {
      /* get vert coord */
      xyz = Mesh_GetVertex( mesh, ii );
      Spherical_XYZ2RTP2D( xyz, rtp );

      area += rtp[0];
      
      /* get vel vec */
      FeVariable_GetValueAtNode( velocityField, ii, vec );
      Spherical_VectorXYZ2RTP( vec, xyz, dim, vrtp );

      // integrate the square
      vrms[0] += (vrtp[0]*vrtp[0]*rtp[0]);
      vrms[1] += (vrtp[1]*vrtp[1]*rtp[0]);
   }

   /* Sum processor integral */
   (void)MPI_Allreduce( vrms, gVrms, 2, MPI_DOUBLE, MPI_SUM, context->communicator );
   (void)MPI_Allreduce( &area, &gArea, 1, MPI_DOUBLE, MPI_SUM, context->communicator );

   // sqrt of volume averaged component
   gVrms[0] = sqrt(gVrms[0]/gArea);
   gVrms[1] = sqrt(gVrms[1]/gArea);

   magVrms = sqrt( gVrms[0]*gVrms[0] + gVrms[1]*gVrms[1] );
   /* Print data to file */
   StgFEM_FrequentOutput_PrintValue( context, magVrms );
   StgFEM_FrequentOutput_PrintValue( context, gVrms[0] );
   StgFEM_FrequentOutput_PrintValue( context, gVrms[1] );

   /* Put Value onto context */
   self->vrms = magVrms;
}

void Underworld_SphericalVrms_PrintHeaderToFile( void* context ) {
   StgFEM_FrequentOutput_PrintString( context, "SphericalVrms" );
   StgFEM_FrequentOutput_PrintString( context, "SphericalVrms_R" );
   StgFEM_FrequentOutput_PrintString( context, "SphericalVrms_T" );
}



