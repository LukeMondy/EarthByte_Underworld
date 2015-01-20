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
#include <Spherical/Spherical.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "CubedSphereVrms.h"

const Type Spherical_CubedSphereVrms_Type = "Spherical_CubedSphereVrms";

void _Spherical_CubedSphereVrms_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
   Spherical_CubedSphereVrms* self = (Spherical_CubedSphereVrms*)component;

   self->context = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"Context", UnderworldContext, True, data );

   self->gaussSwarm = Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"GaussSwarm", Swarm, True, data );

   self->velocityField = Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"VelocityField", FeVariable, True, data );

   Spherical_CubedSphereVrms_PrintHeaderToFile( self->context );

   ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Spherical_CubedSphereVrms_Dump );
}

void _Spherical_CubedSphereVrms_Build( void* component, void* data ) {
   Spherical_CubedSphereVrms* self = (Spherical_CubedSphereVrms*)component;

   assert( self );

   Stg_Component_Build( self->gaussSwarm, data, False );
   Stg_Component_Build( self->velocityField, data, False );
   
   _Codelet_Build( self, data );
}

void _Spherical_CubedSphereVrms_Initialise( void* component, void* data ) {
   Spherical_CubedSphereVrms* self = (Spherical_CubedSphereVrms*)component;

   assert( self );

   Stg_Component_Initialise( self->gaussSwarm, data, False );
   Stg_Component_Initialise( self->velocityField, data, False );
   
   _Codelet_Initialise( self, data );
}

void _Spherical_CubedSphereVrms_Destroy( void* component, void* data ) {
   Spherical_CubedSphereVrms* self = (Spherical_CubedSphereVrms*)component;

   assert( self );

   _Codelet_Destroy( self, data );
   
   Stg_Component_Destroy( self->gaussSwarm, data, False );
   Stg_Component_Destroy( self->velocityField, data, False );
}

void* _Spherical_CubedSphereVrms_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(Spherical_CubedSphereVrms);
   Type                                                      type = Spherical_CubedSphereVrms_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Spherical_CubedSphereVrms_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _Spherical_CubedSphereVrms_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = _Spherical_CubedSphereVrms_Build;
   Stg_Component_InitialiseFunction*                  _initialise = _Spherical_CubedSphereVrms_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Spherical_CubedSphereVrms_Destroy;

   /* 
    * Variables that are set to ZERO are variables that will be set either by the
    * current _New function or another parent _New function further up the hierachy.
    */

    /* default value NON_GLOBAL */
   AllocationType nameAllocationType = NON_GLOBAL;

   return _Codelet_New( CODELET_PASSARGS );
}

Index Spherical_CubedSphereVrms_Register( PluginsManager* pluginsManager ) {
   return PluginsManager_Submit( pluginsManager, Spherical_CubedSphereVrms_Type, (Name)"0", _Spherical_CubedSphereVrms_DefaultNew );
}

/* Integrate Every Step and dump to file */
void Spherical_CubedSphereVrms_Dump( void* _context ) {
   Spherical_CubedSphereVrms* 	self;
   UnderworldContext* 		context 	= (UnderworldContext* ) _context;
   FeMesh*              	mesh;
   FeVariable* 			velocityField	= NULL;
   Swarm* 			is		= NULL;
   IntegrationPoint*		ip		= NULL;
   int 				cell_I, cParticleCount, p_i;
   double 			vec[3], rtp[3], xyz[3], vrtp[3];
   double 			area, detJac, factor, gArea, magVrms;
   double 			vrms[3], gVrms[3]; 
   Dimension_Index		dim 		= context->dim;
   int 				ii, nLocalEls;

   self 	 = (Spherical_CubedSphereVrms*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Spherical_CubedSphereVrms_Type ); assert( self );
   velocityField = self->velocityField; assert( velocityField );
   mesh 	 = velocityField->feMesh; assert( mesh );
   is 		 = self->gaussSwarm;
   nLocalEls 	 = FeMesh_GetElementLocalSize( mesh );
   
   // zero integrals
   vrms[0] = vrms[1] = vrms[2] = area = 0.0;

   for( ii=0; ii<nLocalEls; ii++ ) {
      cell_I = CellLayout_MapElementIdToCellId( is->cellLayout, ii );
      cParticleCount = is->cellParticleCountTbl[ cell_I ];

      for( p_i = 0; p_i < cParticleCount; p_i++ ) {
         ip = (IntegrationPoint*)Swarm_ParticleInCellAt( is, cell_I, p_i );

         /* get vert coord */
         FeMesh_CoordLocalToGlobal( mesh, ii, ip->xi, xyz );
         Spherical_XYZ2regionalSphere( xyz, rtp );

         /* get vel vec */
         FeVariable_InterpolateFromMeshLocalCoord( velocityField, mesh, ii, ip->xi, vec );
         Spherical_VectorXYZ2regionalSphere( vec, xyz, vrtp );

         /* get integration weight and volume */
         detJac = ElementType_JacobianDeterminant( FeMesh_GetElementType( mesh, ii ), mesh, ii, ip->xi, dim );

         factor = detJac*ip->weight;
         area  += factor;

         // integrate the square
         vrms[0] += (vrtp[0]*vrtp[0]*factor);
         vrms[1] += (vrtp[1]*vrtp[1]*factor);
         vrms[2] += (vrtp[2]*vrtp[2]*factor);
      }
   }

   /* Sum processor integral */
   (void)MPI_Allreduce( vrms, gVrms, 3, MPI_DOUBLE, MPI_SUM, context->communicator );
   (void)MPI_Allreduce( &area, &gArea, 1, MPI_DOUBLE, MPI_SUM, context->communicator );

   // sqrt of volume averaged component
   gVrms[0] = sqrt(gVrms[0]/gArea);
   gVrms[1] = sqrt(gVrms[1]/gArea);
   gVrms[2] = sqrt(gVrms[2]/gArea);

   magVrms = sqrt( gVrms[0]*gVrms[0] + gVrms[1]*gVrms[1] + gVrms[2]*gVrms[2] );
   /* Print data to file */
   StgFEM_FrequentOutput_PrintValue( context, magVrms );
   StgFEM_FrequentOutput_PrintValue( context, gVrms[0] );
   StgFEM_FrequentOutput_PrintValue( context, gVrms[1] );
   StgFEM_FrequentOutput_PrintValue( context, gVrms[2] );

   /* Put Value onto context */
   self->vrms = magVrms;
}

void Spherical_CubedSphereVrms_PrintHeaderToFile( void* context ) {
   StgFEM_FrequentOutput_PrintString( context, "CubedSphereVrms" );
   StgFEM_FrequentOutput_PrintString( context, "CubedSphereVrms_R" );
   StgFEM_FrequentOutput_PrintString( context, "CubedSphereVrms_Eta" );
   StgFEM_FrequentOutput_PrintString( context, "CubedSphereVrms_Zeta" );
}
