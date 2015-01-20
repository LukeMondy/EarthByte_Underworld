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
\details CubedSphereNusselt number 
**/
#include <mpi.h>
#include <stdlib.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Spherical/Spherical.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

const Type Spherical_CubedSphereNusselt_Type = "Spherical_CubedSphereNusselt";

typedef struct {
   __Codelet
   Swarm*		gaussSwarm;
   FeMesh*		mesh;
   OperatorFeVariable*	advectiveHeatFluxField;
   OperatorFeVariable*	temperatureTotalDerivField;
   FeVariable*		velocityField;
   FeVariable*		temperatureField;
   FeVariable*		temperatureGradientsField;
   PpcIntegral*		volAvgVD;
   PpcIntegral*		volAvgW;
   PpcIntegral*		volAvgT;
   PpcIntegral*		vol;
   PpcManager*		mgr;
   double		Ra;
} Spherical_CubedSphereNusselt;

Index Spherical_CubedSphereNusselt_Register( PluginsManager* pluginsManager );
void Spherical_MaxVel_Output( UnderworldContext* context ) ;
void Spherical_CubedSphereNusselt_Output( UnderworldContext* context ) ;
void Spherical_Work_Output( UnderworldContext* context ) ;

void _Spherical_CubedSphereNusselt_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
   Spherical_CubedSphereNusselt*	self 	= (Spherical_CubedSphereNusselt*)component;
   UnderworldContext*			context;
   FieldVariable_Register*		fV_Register;
   FeVariable*				temperatureGradientsField;
   FeVariable*				velocityField;
   FeVariable*				temperatureField;

   self->context = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"Context", UnderworldContext, True, data );

   self->gaussSwarm = Stg_ComponentFactory_PluginConstructByKey(
      cf, self, (Dictionary_Entry_Key)"GaussSwarm", Swarm, True, data );

   /* The PpcManager */
   self->mgr = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
   if( !self->mgr  )
      self->mgr = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, True, data  );

   // initialise as unfound ppc
   self->volAvgVD = self->volAvgW = self->volAvgT = self->vol = NULL;

   self->volAvgVD = Stg_ComponentFactory_PluginConstructByKey( cf, self, "volume_averaged_viscous_dissipation",PpcIntegral, False, data);
   self->volAvgW = Stg_ComponentFactory_PluginConstructByKey( cf, self, "volume_averaged_work_done",PpcIntegral, False, data);
   self->volAvgT = Stg_ComponentFactory_PluginConstructByKey( cf, self, "volume_averaged_temperature",PpcIntegral, True, data);
   self->vol = Stg_ComponentFactory_PluginConstructByKey( cf, self, "volume",PpcIntegral, True, data);

   self->Ra = Stg_ComponentFactory_PluginGetDouble( cf, self, (Dictionary_Entry_Key)"Ra", -2 );
   assert( self->Ra > -1 );

   StgFEM_FrequentOutput_PrintString( self->context, "<T'>" );
   StgFEM_FrequentOutput_PrintString( self->context, "NuU_av" );
   StgFEM_FrequentOutput_PrintString( self->context, "NuB_av" );

   StgFEM_FrequentOutput_PrintString( self->context, "out_max_VelMag" );
   StgFEM_FrequentOutput_PrintString( self->context, "in_max_VelMag" );

   // if ppcs found then add to FrequentOutput file
   if( self->vol ) {
      if( self->volAvgVD )
         StgFEM_FrequentOutput_PrintString( self->context, "volAvgVD" );
      if( self->volAvgW )
         StgFEM_FrequentOutput_PrintString( self->context, "volAvgW" );
   }
   
   context = (UnderworldContext*)self->context;
   fV_Register = context->fieldVariable_Register;
   /* Create Some FeVariables to calculate nusselt number */
   temperatureField = self->temperatureField = (FeVariable*)FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
   velocityField = self->velocityField = (FeVariable*)FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
   temperatureGradientsField = self->temperatureGradientsField = (FeVariable*)FieldVariable_Register_GetByName( fV_Register, "TemperatureGradientsField" );
   
   self->mesh = ((FeVariable*)temperatureField)->feMesh;
   self->advectiveHeatFluxField = OperatorFeVariable_NewBinary(  
      "AdvectiveHeatFluxField", (DomainContext*)context, temperatureField, velocityField, "VectorScale" );

   self->temperatureTotalDerivField = OperatorFeVariable_NewBinary(  
      "TemperatureTotalDerivField", (DomainContext*)context, self->advectiveHeatFluxField, temperatureGradientsField, 
      "Subtraction" );
   
   /* Add functions to entry points */
   ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Spherical_CubedSphereNusselt_Output );
   ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Spherical_MaxVel_Output );

   // if the definition for the vd and adiabatic work exist add another hook
   if( self->vol && self->volAvgVD && self->volAvgW )
      ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Spherical_Work_Output );
}

void _Spherical_CubedSphereNusselt_Build( void* component, void* data ) {
   Spherical_CubedSphereNusselt*   self = (Spherical_CubedSphereNusselt*)component;

   Stg_Component_Build( self->gaussSwarm, data, False );
   Stg_Component_Build( self->advectiveHeatFluxField, data, False );
   Stg_Component_Build( self->temperatureTotalDerivField, data, False );
   
   _Codelet_Build( component, data );
}

void _Spherical_CubedSphereNusselt_Initialise( void* component, void* data ) {
   Spherical_CubedSphereNusselt*   self = (Spherical_CubedSphereNusselt*)component;

   Stg_Component_Initialise( self->gaussSwarm, data, False );
   Stg_Component_Initialise( self->advectiveHeatFluxField, data, False );
   Stg_Component_Initialise( self->temperatureTotalDerivField, data, False );
   
   _Codelet_Initialise( component, data );
}

void _Spherical_CubedSphereNusselt_Destroy( void* component, void* data ) {
   Spherical_CubedSphereNusselt*   self = (Spherical_CubedSphereNusselt*)component;

   Stg_Component_Destroy( self->gaussSwarm, data, False );
   Stg_Component_Destroy( self->advectiveHeatFluxField, data, False );
   Stg_Component_Destroy( self->temperatureTotalDerivField, data, False );
   
   _Codelet_Destroy( component, data );
}

void* _Spherical_CubedSphereNusselt_DefaultNew( Name name ) {
   /* Variables set in this function */
   SizeT                                             _sizeOfSelf = sizeof(Spherical_CubedSphereNusselt);
   Type                                                     type = Spherical_CubedSphereNusselt_Type;
   Stg_Class_DeleteFunction*                             _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                               _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                 _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction* _defaultConstructor = _Spherical_CubedSphereNusselt_DefaultNew;
   Stg_Component_ConstructFunction*                   _construct = _Spherical_CubedSphereNusselt_AssignFromXML;
   Stg_Component_BuildFunction*                           _build = _Spherical_CubedSphereNusselt_Build;
   Stg_Component_InitialiseFunction*                 _initialise = _Spherical_CubedSphereNusselt_Initialise;
   Stg_Component_ExecuteFunction*                       _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                       _destroy = _Spherical_CubedSphereNusselt_Destroy;

   /* 
    * Variables that are set to ZERO are variables that will be set either by the
    * current _New function or another parent _New function further up the hierachy.
    */
   AllocationType  nameAllocationType = NON_GLOBAL;

   return _Codelet_New( CODELET_PASSARGS );
}

Index Spherical_CubedSphereNusselt_Register( PluginsManager* pluginsManager ) {
   return PluginsManager_Submit( pluginsManager, Spherical_CubedSphereNusselt_Type, (Name)"0", _Spherical_CubedSphereNusselt_DefaultNew );
}

void Spherical_Work_Output( UnderworldContext* context ) {
   Spherical_CubedSphereNusselt* self=(Spherical_CubedSphereNusselt*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Spherical_CubedSphereNusselt_Type );
   double vol, volAvgVD, volAvgW;

   // get energy metrics
   vol = PpcIntegral_Integrate( self->vol );
   volAvgVD = PpcIntegral_Integrate( self->volAvgVD );
   volAvgW = PpcIntegral_Integrate( self->volAvgW );

   // make volume averaged values
   volAvgVD = volAvgVD/vol;
   volAvgW = volAvgW/vol;

   //print
   StgFEM_FrequentOutput_PrintValue( context, volAvgVD ); 
   StgFEM_FrequentOutput_PrintValue( context, volAvgW );
}
void Spherical_MaxVel_Output( UnderworldContext* context ) {
   Spherical_CubedSphereNusselt* 	self		= (Spherical_CubedSphereNusselt*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Spherical_CubedSphereNusselt_Type );
   FeVariable*				velocityField   = self->velocityField;
   FeMesh* 				mesh 		= velocityField->feMesh;
   Grid*				grid		= NULL;
   unsigned*				sizes		= NULL;
   int 					vertId, dId;
   unsigned 				dT_i, dT_j, nDomainSize, ijk[3];
   double 				vel[3], velMax[2], gVelMax[2], velMag;

   //set 'em big and negative
   velMax[0] = velMax[1] = -1*HUGE_VAL;

   // get vert grid
   RegularMeshUtils_ErrorCheckAndGetDetails( (Mesh*)mesh, MT_VERTEX, &nDomainSize, &grid );
   sizes = grid->sizes;

   // go around 
   for( dT_i = 0; dT_i < sizes[1]; dT_i++ ) {
      for( dT_j = 0; dT_j < sizes[2]; dT_j++ ) {
         // find inner vertex
         ijk[0] = 0;
         // angular discretisation
         ijk[1] = dT_i;
         ijk[2] = dT_j;
         vertId = Grid_Project( grid, ijk );

         // if the node is local
         if( Mesh_GlobalToDomain( mesh, MT_VERTEX, vertId, &dId ) == True ) {
            FeVariable_GetValueAtNode( velocityField, dId, vel );
            velMag = sqrt( vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2] );
            if( velMag > velMax[0] ) velMax[0] = velMag;
         }

         // find outer vertex
         ijk[0] = grid->sizes[0]-1;
         vertId = Grid_Project( grid, ijk );

         // if the node is local
         if( Mesh_GlobalToDomain( mesh, MT_VERTEX, vertId, &dId ) == True ) {
            FeVariable_GetValueAtNode( velocityField, dId, vel );
            velMag = sqrt( vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2] );
            if( velMag > velMax[1] ) velMax[1] = velMag;
         }
      }
   }

   (void)MPI_Allreduce( velMax, gVelMax, 2, MPI_DOUBLE, MPI_MAX, context->communicator );

   StgFEM_FrequentOutput_PrintValue( context, gVelMax[1] ); // print outer velocity max
   StgFEM_FrequentOutput_PrintValue( context, gVelMax[0] ); // print inner velocity max
}
      
void Spherical_CubedSphereNusselt_Output( UnderworldContext* context ) {
   Spherical_CubedSphereNusselt* 	self  		= NULL;
   Swarm*				swarm		= NULL;
   FeMesh* 				mesh		= NULL;
   ElementType* 			elementType	= NULL;
   Grid*				grid		= NULL;
   double 				avgT, vol;

   self = (Spherical_CubedSphereNusselt*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Spherical_CubedSphereNusselt_Type );

   FeVariable 	*temperatureGradientsField 	= self->temperatureGradientsField;
   FeVariable 	*temperatureField 		= self->temperatureField;
   FeVariable 	*velocityField 			= self->velocityField;

   swarm = self->gaussSwarm;
   mesh = (FeMesh*)self->mesh;

   assert( self );
   assert( mesh );
   assert( swarm );

   /*
      for each boundary element at the top
         integrate dT/dr &
         integrate T*v_r

      ASSUMPTIONS:
      * uses grid data structure for r,theta mesh (cartesian topology so i can)
   */

   unsigned 		nEls, e_i, cell_i, nPoints, p_i, ijk[3];
   IntegrationPoint*	particle;
   double 		value[3], xyz[3], rtp[3], detJac, dT_dr, factor, vel[3], T, vel_rtp[3];
   double 		gVolume[2], volume[2];
   double 		J_Nu[2], gJ_Nu[2];
   double		rMin			= ((RSGenerator*)mesh->generator)->crdMin[0];
   double		rMax			= ((RSGenerator*)mesh->generator)->crdMin[1];
   double		dxdr[3], r;

   elementType = FeMesh_GetElementType( mesh, 0 ); // assuming all element are the same as el 0

   memset( J_Nu, 0, 3*sizeof(double) );
   memset( volume, 0, 3*sizeof(double) );

   RegularMeshUtils_ErrorCheckAndGetDetails( (Mesh*)mesh, MT_VOLUME, &nEls, &grid );

   for( e_i = 0; e_i < nEls; e_i++ ) {
      // use cartesian grid data structure - can improve later on to make more general
      RegularMeshUtils_Element_1DTo3D( mesh, Mesh_DomainToGlobal( (Mesh*)mesh, MT_FACE, e_i ), ijk );

      // if element is on the outer radius boundary
      if( ijk[0] == grid->sizes[0] - 1 ) {
         // get integration number of integration points in cell
         cell_i = CellLayout_MapElementIdToCellId( swarm->cellLayout, e_i );
         nPoints = swarm->cellParticleCountTbl[ cell_i ];

         for( p_i = 0; p_i < nPoints; p_i++ ) {
            // get integration particle
            particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_i, p_i );

            // get temperatureDeriv and xyz
            FeVariable_InterpolateFromMeshLocalCoord( temperatureField, (FeMesh*)mesh, e_i, particle->xi, &T);
            FeVariable_InterpolateFromMeshLocalCoord( velocityField, (FeMesh*)mesh, e_i, particle->xi, vel);
            FeVariable_InterpolateFromMeshLocalCoord( temperatureGradientsField, (FeMesh*)mesh, e_i, particle->xi, value);
            FeMesh_CoordLocalToGlobal( mesh, e_i, particle->xi, xyz );

            Spherical_XYZ2regionalSphere( xyz, rtp );
            Spherical_VectorXYZ2regionalSphere( vel, xyz, vel_rtp );
            detJac = ElementType_JacobianDeterminant( elementType, (FeMesh*)mesh, e_i, particle->xi, 3 );

            r = sqrt( rtp[0]*rtp[0] + rtp[1]*rtp[1] + rtp[2]*rtp[2] );
            dxdr[0] = r/xyz[0];
            dxdr[1] = r/xyz[1];
            dxdr[2] = r/xyz[2];
            // calc dT_dr = dT_dx * dx_dr + dT_dy * dy_dr + dT_dz * dz_dr
            dT_dr = value[0]*dxdr[0] + value[1]*dxdr[1] + value[2]*dxdr[2];

            // add to element integral
            J_Nu[0] += particle->weight * detJac * (dT_dr - T*vel_rtp[0]);
            volume[0] += detJac * particle->weight;
         }
      }
      // if element is on the inner radius boundary
      if( ijk[0] == 0 ) {
         // get integration number of integration points in cell
         cell_i = CellLayout_MapElementIdToCellId( swarm->cellLayout, e_i );
         nPoints = swarm->cellParticleCountTbl[ cell_i ];

         for( p_i = 0; p_i < nPoints; p_i++ ) {
            // get integration particle
            particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_i, p_i );

            // get temperatureDeriv and xyz
            FeVariable_InterpolateFromMeshLocalCoord( temperatureField, (FeMesh*)mesh, e_i, particle->xi, &T);
            FeVariable_InterpolateFromMeshLocalCoord( velocityField, (FeMesh*)mesh, e_i, particle->xi, vel);
            FeVariable_InterpolateFromMeshLocalCoord( temperatureGradientsField, (FeMesh*)mesh, e_i, particle->xi, value);
            FeMesh_CoordLocalToGlobal( mesh, e_i, particle->xi, xyz );

            Spherical_XYZ2regionalSphere( xyz, rtp );
            Spherical_VectorXYZ2regionalSphere( vel, xyz, vel_rtp );
            detJac = ElementType_JacobianDeterminant( elementType, (FeMesh*)mesh, e_i, particle->xi, 3 );

            r = sqrt( rtp[0]*rtp[0] + rtp[1]*rtp[1] + rtp[2]*rtp[2] );
            dxdr[0] = r/xyz[0];
            dxdr[1] = r/xyz[1];
            dxdr[2] = r/xyz[2];
            // calc dT_dr = dT_dx * dx_dr + dT_dy * dy_dr + dT_dz * dz_dr
            dT_dr = value[0]*dxdr[0] + value[1]*dxdr[1] + value[2]*dxdr[2];

            J_Nu[1] += particle->weight * detJac * (dT_dr - T*vel_rtp[0]);
            volume[1] += detJac * particle->weight;
         }
      }
   }

   /* Sum of procs integral */
   (void)MPI_Allreduce( J_Nu, gJ_Nu, 2, MPI_DOUBLE, MPI_SUM, context->communicator );
   (void)MPI_Allreduce( volume, gVolume, 2, MPI_DOUBLE, MPI_SUM, context->communicator );

   // to get horizontally averaged quantities we divide by volume
   gJ_Nu[0] /= gVolume[0];
   gJ_Nu[1] /= gVolume[1];

   // normalise CubedSphereNusselt upper condition - this is for scaling against published results
   // 1.22 and 2.22 are the inner and outer radii for those results
   factor = rMax * log(0.55);
   gJ_Nu[0] = factor *  gJ_Nu[0];

   // normalise CubedSphereNusselt lower condition
   factor = rMin * log(0.55);
   gJ_Nu[1] = factor * gJ_Nu[1];

   avgT = PpcIntegral_Integrate( self->volAvgT );
   vol = PpcIntegral_Integrate( self->vol );

   StgFEM_FrequentOutput_PrintValue( context, avgT/vol );
   StgFEM_FrequentOutput_PrintValue( context, gJ_Nu[0] );
   StgFEM_FrequentOutput_PrintValue( context, gJ_Nu[1] );
}
