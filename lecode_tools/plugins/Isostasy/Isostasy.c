/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
** Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
** Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
** AuScope - http://www.auscope.org
** Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
** Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
** Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
** Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
** Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
** David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
** Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
** Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
** Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
** Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
** Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
** Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
** Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
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

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <stddef.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include "Isostasy.h"


const Type lecode_tools_Isostasy_Type = "lecode_tools_Isostasy";
lecode_tools_Isostasy *lecode_tools_Isostasy_self;


void lecode_tools_Isostasy_AverageSurfaces(lecode_tools_Isostasy *self,
      double** _avg_sep, double** _avg_height)
{
   FeMesh *mesh;
   Grid *grid, *elGrid;
   double *local_height, *global_height;
   double *local_top_vy, *global_top_vy;
   double *local_bot_vy, *global_bot_vy;
   int *local_top_cnt, *global_top_cnt;
   int *local_bot_cnt, *global_bot_cnt;
   int param[3], elParam[3], num_nodes, n;
   double vel[3];
   IArray *inc;
   int nDims, arrayPos, arraySize;
   int ii, jj;

   mesh = self->mesh;
   nDims = Mesh_GetDimSize( mesh );
   grid = *Mesh_GetExtension(mesh, Grid**, mesh->vertGridId );
   elGrid = *Mesh_GetExtension(mesh, Grid**, mesh->elGridId );
   num_nodes = FeMesh_GetNodeLocalSize(mesh);
   inc = IArray_New();

   /* mem alloc from bottom surface */
   arraySize=0; /*to prevent warnings*/
   if ( nDims == 2 ) arraySize = elGrid->sizes[0];
   else if ( nDims == 3 ) arraySize = elGrid->sizes[0]*elGrid->sizes[self->zontalAxis];
   else assert(0);

   local_top_vy = (double*)malloc( arraySize*sizeof(double) );
   memset( local_top_vy, 0, arraySize*sizeof(double) );
   local_bot_vy = (double*)malloc( arraySize*sizeof(double) );
   memset( local_bot_vy, 0, arraySize*sizeof(double) );
   local_height = (double*)malloc( arraySize*sizeof(double) );
   memset( local_height, 0, arraySize*sizeof(double) );
   local_top_cnt = (int*)malloc( arraySize*sizeof(int) );
   memset( local_top_cnt, 0, arraySize*sizeof(int) );
   local_bot_cnt = (int*)malloc( arraySize*sizeof(int) );
   memset( local_bot_cnt, 0, arraySize*sizeof(int) );

   for (ii = 0; ii < num_nodes; ii++)
   {

      FeMesh_GetNodeElements( mesh, ii, inc );

      n = FeMesh_NodeDomainToGlobal(mesh, ii);
      Grid_Lift(grid, n, param);

      if ((self->surfaceIdx != -1 && param[self->vertAxis] == self->surfaceIdx) ||
            (self->surfaceIdx == -1 && param[self->vertAxis] == grid->sizes[self->vertAxis] - 1))
      {
         FeVariable_GetValueAtNode(self->vel_field, ii, vel);

         if ( self->avg )
         {
            local_top_vy[0] += vel[self->vertAxis];
            local_height[0] += Mesh_GetVertex( mesh, ii )[self->vertAxis];
            local_top_cnt[0]++;
         }
         else
         {
            for (jj = 0; jj < inc->size; jj++ )
            {
               Grid_Lift( elGrid, FeMesh_ElementDomainToGlobal( self->mesh, inc->ptr[jj] ), elParam );

               /* Make sure element is below surface. */
               if ( self->surfaceIdx != -1 && elParam[self->vertAxis] >= self->surfaceIdx )
                  continue;

               arrayPos = elParam[0];
               if ( nDims == 3 ) arrayPos += elParam[self->zontalAxis]*elGrid->sizes[0];

               local_top_vy[arrayPos] += vel[self->vertAxis];
               local_height[arrayPos] += Mesh_GetVertex( mesh, ii )[self->vertAxis];
               local_top_cnt[arrayPos]++;
            }
         }
      }

      if (param[self->vertAxis] == 0 )
      {
         FeVariable_GetValueAtNode(self->vel_field, ii, vel);

         if ( self->avg )
         {
            local_bot_vy[0] += vel[self->vertAxis];
            local_height[0] -= Mesh_GetVertex( mesh, ii )[self->vertAxis];
            local_bot_cnt[0]++;
         }
         else
         {
            for (jj = 0; jj < inc->size; jj++ )
            {
               Grid_Lift( elGrid, FeMesh_ElementDomainToGlobal( self->mesh, inc->ptr[jj] ), elParam );

               /* Make sure element is below surface. */
               if ( self->surfaceIdx != -1 && elParam[self->vertAxis] >= self->surfaceIdx )
                  continue;

               arrayPos = elParam[0];
               if ( nDims == 3 ) arrayPos += elParam[self->zontalAxis]*elGrid->sizes[0];

               local_bot_vy[arrayPos] += vel[self->vertAxis];
               local_height[arrayPos] -= Mesh_GetVertex( mesh, ii )[self->vertAxis];
               local_bot_cnt[arrayPos]++;
            }
         }
      }

   }

   global_top_vy = (double*)malloc( arraySize*sizeof(double) );
   global_bot_vy = (double*)malloc( arraySize*sizeof(double) );
   global_height = (double*)malloc( arraySize*sizeof(double) );
   global_top_cnt = (int*)malloc( arraySize*sizeof(int) );
   global_bot_cnt = (int*)malloc( arraySize*sizeof(int) );

   MPI_Allreduce(local_height, global_height, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_top_vy, global_top_vy, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_bot_vy, global_bot_vy, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_top_cnt, global_top_cnt, arraySize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_bot_cnt, global_bot_cnt, arraySize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   free( local_height );
   free( local_top_vy );
   free( local_bot_vy );
   free( local_top_cnt );
   free( local_bot_cnt );
   Stg_Class_Delete( inc );

   *_avg_sep = (double*)malloc( arraySize*sizeof(double) );
   *_avg_height = (double*)malloc( arraySize*sizeof(double) );
   for ( ii = 0; ii < arraySize; ii++ )
   {
      (*_avg_sep)[ii] = global_top_vy[ii]/(double)(global_top_cnt[ii]) - global_bot_vy[ii]/(double)(global_bot_cnt[ii]);
      (*_avg_height)[ii] = global_height[ii]/(double)(global_top_cnt[ii]);
      if (self->avg)
         break;
   }

   free( global_height );
   free( global_top_vy );
   free( global_bot_vy );
   free( global_top_cnt );
   free( global_bot_cnt );
}


void lecode_tools_Isostasy_AverageBody(lecode_tools_Isostasy *self,
                                       double** _avg_density, double** _rho_zero_density,
                                       double** _phi)
{
   FeMesh *mesh;
   ElementType *el_type;
   IntegrationPointsSwarm *swarm;
   double *local_density, *global_density;
   double *local_vol, *global_vol;
   double *local_rho_zero_vol, *global_rho_zero_vol;
   double *local_rho_zero_density, *global_rho_zero_density;
   double *local_phi, *global_phi, temp, tempDot;
   int cell, num_particles, num_dims, num_els;
   IntegrationPoint *particle;
   double jac_det;
   double density, alpha, densityFinal;
   Material *mat;
   Bool oneToMany;
   Grid* elGrid;
   int elInds[3], arraySize, arrayPos;
   int ii, jj;

   mesh = self->mesh;
   elGrid = *Mesh_GetExtension( mesh, Grid**,  mesh->elGridId );
   num_dims = Mesh_GetDimSize(mesh);
   swarm = self->swarm;
   num_els = FeMesh_GetElementLocalSize(mesh);

   arraySize=0;
   if ( num_dims == 2 ) arraySize = elGrid->sizes[0];
   else if ( num_dims == 3 ) arraySize = elGrid->sizes[0]*elGrid->sizes[self->zontalAxis];
   else assert(0);

   /* Allocate for the column values. */
   local_vol = (double*)malloc( arraySize*sizeof(double) );
   memset( local_vol, 0, arraySize*sizeof(double) );
   local_density = (double*)malloc( arraySize*sizeof(double) );
   memset( local_density, 0, arraySize*sizeof(double) );
   local_rho_zero_vol = (double*)malloc( arraySize*sizeof(double) );
   memset( local_rho_zero_vol, 0, arraySize*sizeof(double) );
   local_rho_zero_density = (double*)malloc( arraySize*sizeof(double) );
   memset( local_rho_zero_density, 0, arraySize*sizeof(double) );
   local_phi = (double*)malloc( arraySize*sizeof(double) );
   memset( local_phi, 0, arraySize*sizeof(double) );

   /* Initialise temperature. */
   temp = 0.0;

   oneToMany = Stg_Class_IsInstance(swarm->mapper, OneToManyMapper_Type);

   for (ii = 0; ii < num_els; ii++)
   {

      /* Make sure the element is beneath the surface. */
      Grid_Lift( elGrid, FeMesh_ElementDomainToGlobal( mesh, ii ), elInds );
      if ( self->surfaceIdx != -1 && elInds[self->vertAxis] >= self->surfaceIdx )
         continue;

      el_type = FeMesh_GetElementType(mesh, ii);
      cell = CellLayout_MapElementIdToCellId(swarm->cellLayout, ii);
      num_particles = swarm->cellParticleCountTbl[cell];

      for (jj = 0; jj < num_particles; jj++)
      {

         particle = (IntegrationPoint*)Swarm_ParticleInCellAt(swarm, cell, jj);
         jac_det = ElementType_JacobianDeterminant(el_type, mesh, ii, particle->xi, num_dims);

         if(!self->ppcManager){
            density = IntegrationPointMapper_GetDoubleFromMaterial(
                         swarm->mapper, particle, self->buoyancy->materialExtHandle,
                         offsetof(BuoyancyForceTerm_MaterialExt, density) );
            alpha = IntegrationPointMapper_GetDoubleFromMaterial(
                       swarm->mapper, particle, self->buoyancy->materialExtHandle,
                       offsetof(BuoyancyForceTerm_MaterialExt, alpha) );

            if (self->tempField)
            {
               FeVariable_InterpolateFromMeshLocalCoord(self->tempField, self->tempField->feMesh,
                     ii, particle->xi, &temp);
               FeVariable_InterpolateFromMeshLocalCoord(self->tempDotField, self->tempDotField->feMesh,
                     ii, particle->xi, &tempDot);
            }

            densityFinal = density*(1.0 - alpha*temp);

         } else {
            int err;
            /* Density */
            err = PpcManager_Get( self->ppcManager, ii, particle, self->densityID, &densityFinal );
            assert(!err);
         }

         arrayPos = elInds[0];
         if ( num_dims == 3 ) arrayPos += elInds[self->zontalAxis]*elGrid->sizes[0];

         local_vol[arrayPos] += particle->weight*jac_det;
         local_density[arrayPos] += particle->weight*jac_det*densityFinal;

         if (!oneToMany)
         {
            mat = IntegrationPointsSwarm_GetMaterialOn(swarm, particle);
            if (mat->index == self->rho_zero_mat->index)
            {
               local_rho_zero_vol[arrayPos] += particle->weight*jac_det;
               local_rho_zero_density[arrayPos] += particle->weight*jac_det*densityFinal;
            }
         }
         else
         {
            OneToManyRef *ref;
            int cnt;
            int kk;

            ref = OneToManyMapper_GetMaterialRef(swarm->mapper, particle);
            cnt = 0;
            for (kk = 0; kk < ref->numParticles; kk++)
            {
               mat = MaterialPointsSwarm_GetMaterialAt(((OneToManyMapper*)swarm->mapper)->materialSwarm, ref->particleInds[kk]);
               if (mat->index == self->rho_zero_mat->index)
                  cnt++;
            }

            if (2*cnt > ref->numParticles)
            {
               local_rho_zero_vol[arrayPos] += particle->weight*jac_det;
               local_rho_zero_density[arrayPos] += particle->weight*jac_det*densityFinal;
            }
         }

         if (_phi)
         {
            local_phi[arrayPos] += particle->weight*jac_det*(-density*alpha*tempDot);
         }
      }
   }

   /* Allocate for the global column values. */
   global_vol = (double*)malloc( arraySize*sizeof(double) );
   global_density = (double*)malloc( arraySize*sizeof(double) );
   global_rho_zero_vol = (double*)malloc( arraySize*sizeof(double) );
   global_rho_zero_density = (double*)malloc( arraySize*sizeof(double) );
   global_phi = (double*)malloc( arraySize*sizeof(double) );

   MPI_Allreduce(local_vol, global_vol, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_density, global_density, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_rho_zero_vol, global_rho_zero_vol, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(local_rho_zero_density, global_rho_zero_density, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (_phi)
      MPI_Allreduce(local_phi, global_phi, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   free( local_vol );
   free( local_density );
   free( local_rho_zero_vol );
   free( local_rho_zero_density );
   free( local_phi );

   if ( self->avg )
   {
      for ( ii = 1; ii < arraySize; ii++ )
      {
         global_vol[0] += global_vol[ii];
         global_density[0] += global_density[ii];
         global_rho_zero_vol[0] += global_rho_zero_vol[ii];
         global_rho_zero_density[0] += global_rho_zero_density[ii];
         if ( _phi )
            global_phi[0] += global_phi[ii];
      }
   }

   /* Calculate results. */
   *_avg_density = (double*)malloc( arraySize*sizeof(double) );
   *_rho_zero_density = (double*)malloc( arraySize*sizeof(double) );
   if (_phi)
      *_phi = (double*)malloc( arraySize*sizeof(double) );
   for ( ii = 0; ii < arraySize; ii++ )
   {
      (*_avg_density)[ii] = (global_vol[ii] > 1e-7) ? global_density[ii]/global_vol[ii] : 0.0;
      (*_rho_zero_density)[ii] = (global_rho_zero_vol[ii] > 1e-7) ? global_rho_zero_density[ii]/global_rho_zero_vol[ii] : 0.0;
      if (_phi)
         (*_phi)[ii] = (global_vol[ii] > 1e-7) ? global_phi[ii]/global_vol[ii] : 0.0;
      if ( self->avg )
         break;
   }

   /*
       printf("Global mean density: %g\n", (*_avg_density)[0]);
       printf("Global mean rho_0 density: %g\n", (*_rho_zero_density)[0]);
       printf("Global phi/vol: %g\n", (*_phi)[0]);
   */

   free( global_vol );
   free( global_density );
   free( global_rho_zero_vol );
   free( global_rho_zero_density );
   free( global_phi );
}

void lecode_tools_Isostasy_SolveThermal(lecode_tools_Isostasy *self)
{
   _SystemLinearEquations_Execute(self->thermalSLE, self->ctx);
}


void lecode_tools_Isostasy_CalcBasalFlow(lecode_tools_Isostasy *self)
{
   int rank;
   Grid *nodeGrid, *elGrid;
   double *avg_sep, *avg_height;
   double *avg_density, *rho_zero_density;
   double *phi, **phi_ptr;
   double *node_avg_sep, *node_avg_height;
   double *node_avg_density, *node_rho_zero_density, *node_phi;
   double *global_node_avg_sep, *global_node_avg_height;
   double *global_node_avg_density, *global_node_rho_zero_density, *global_node_phi;
   double avgFlow, tmp;
   int param[3], nodeIdx, arraySize, arrayPos, elPos;
   int nDims, iterSizes[2];
   IArray *inc;
   int ii, jj, kk;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   nDims = Mesh_GetDimSize( self->mesh );
   nodeGrid = *Mesh_GetExtension( self->mesh, Grid**,  self->mesh->vertGridId );
   elGrid = *Mesh_GetExtension( self->mesh, Grid**,  self->mesh->elGridId );
   inc = IArray_New();

   if (self->thermalSLE)
   {
      lecode_tools_Isostasy_SolveThermal(self);
      phi_ptr = &phi;
   }
   else
      phi_ptr = NULL;

   lecode_tools_Isostasy_AverageSurfaces(self, &avg_sep, &avg_height);
   lecode_tools_Isostasy_AverageBody(self, &avg_density, &rho_zero_density, phi_ptr);
   // avg_density seems to be the avg rho across the basal node columns.
   // rho_zero_density seems to be the avg density of the rho_zero material.

   if ( self->avg )
   {

      self->flow[0] = -1.0*avg_density[0]*avg_sep[0]/rho_zero_density[0];
      if (phi_ptr)
         self->flow[0] -= avg_height[0]*phi[0]/rho_zero_density[0];

      avgFlow = self->flow[0];

   }
   else
   {

      // nodeGrid->sizes is the dimensions of the mesh in x (0), y (1), and z(2)
      arraySize=0;
      if ( nDims == 2 ) arraySize = nodeGrid->sizes[0];
      else if ( nDims == 3 ) arraySize = nodeGrid->sizes[0]*nodeGrid->sizes[self->zontalAxis];
      else assert(0);
      // the arraysize is the base of the model, in nodes.

      node_avg_density = (double*)malloc( arraySize*sizeof(double) );
      memset( node_avg_density, 0, arraySize*sizeof(double) );
      node_rho_zero_density = (double*)malloc( arraySize*sizeof(double) );
      memset( node_rho_zero_density, 0, arraySize*sizeof(double) );
      node_avg_sep = (double*)malloc( arraySize*sizeof(double) );
      memset( node_avg_sep, 0, arraySize*sizeof(double) );
      node_avg_height = (double*)malloc( arraySize*sizeof(double) );
      memset( node_avg_height, 0, arraySize*sizeof(double) );
      if ( phi_ptr )
      {
         node_phi = (double*)malloc( arraySize*sizeof(double) );
         memset( node_phi, 0, arraySize*sizeof(double) );
      }
     
      param[0] = param[1] = param[2] = 0;
      iterSizes[0] = nodeGrid->sizes[0];
      iterSizes[1] = (nDims == 3) ? nodeGrid->sizes[self->zontalAxis] : 1;
      for ( ii = 0; ii < iterSizes[1]; ii++ )
      {
         // Going into the 'Y' or 'Z' of the basal surface. If 2D, then this loop only goes once.
         for ( kk = 0; kk < iterSizes[0]; kk++)
         {

            param[self->zontalAxis] = ii;
            param[0] = kk;
            arrayPos = param[0] + param[self->zontalAxis]*nodeGrid->sizes[0];

            nodeIdx = Grid_Project( nodeGrid, param );
            if ( !FeMesh_NodeGlobalToDomain( self->mesh, nodeIdx, &nodeIdx ) ||
                  nodeIdx >= FeMesh_GetNodeLocalSize( self->mesh ) )
            {
               continue;
            }

            FeMesh_GetNodeElements( self->mesh, nodeIdx, inc );
            for ( jj = 0; jj < inc->size; jj++ )
            {
               Grid_Lift( elGrid, FeMesh_ElementDomainToGlobal( self->mesh, inc->ptr[jj] ), param );

               elPos = param[0];
               if ( nDims == 3 ) elPos += param[self->zontalAxis]*elGrid->sizes[0];

               node_avg_density[arrayPos] += avg_density[elPos];
               node_rho_zero_density[arrayPos] += rho_zero_density[elPos];
               node_avg_sep[arrayPos] += avg_sep[elPos];
               node_avg_height[arrayPos] += avg_height[elPos];
               if ( phi_ptr )
                  node_phi[arrayPos] += phi[elPos];
            }

            node_avg_density[arrayPos] /= (double)inc->size;
            node_rho_zero_density[arrayPos] /= (double)inc->size;
            node_avg_sep[arrayPos] /= (double)inc->size;
            node_avg_height[arrayPos] /= (double)inc->size;
            if ( phi_ptr )
               node_phi[arrayPos] /= (double)inc->size;

         }

      }

      global_node_avg_density = (double*)malloc( arraySize*sizeof(double) );
      memset( global_node_avg_density, 0, arraySize*sizeof(double) );
      global_node_rho_zero_density = (double*)malloc( arraySize*sizeof(double) );
      memset( global_node_rho_zero_density, 0, arraySize*sizeof(double) );
      global_node_avg_sep = (double*)malloc( arraySize*sizeof(double) );
      memset( global_node_avg_sep, 0, arraySize*sizeof(double) );
      global_node_avg_height = (double*)malloc( arraySize*sizeof(double) );
      memset( global_node_avg_height, 0, arraySize*sizeof(double) );
      if ( phi_ptr )
      {
         global_node_phi = (double*)malloc( arraySize*sizeof(double) );
         memset( global_node_phi, 0, arraySize*sizeof(double) );
      }

      MPI_Allreduce(node_avg_density, global_node_avg_density, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(node_rho_zero_density, global_node_rho_zero_density, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(node_avg_sep, global_node_avg_sep, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(node_avg_height, global_node_avg_height, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if ( phi_ptr )
         MPI_Allreduce(node_phi, global_node_phi, arraySize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      free( avg_sep );
      free( avg_height );
      free( avg_density );
      free( rho_zero_density );
      if ( phi_ptr )
         free( phi );
      Stg_Class_Delete( inc );

      avgFlow = 0.0;

      for ( ii = 0; ii < arraySize; ii++ )
      {

         /*
           param[0] = ii;
           nodeIdx = Grid_Project( nodeGrid, param );
           if( !FeMesh_NodeGlobalToDomain( self->mesh, nodeIdx, &nodeIdx ) ||
           nodeIdx >= FeMesh_GetNodeLocalSize( self->mesh ) )
           {
           continue;
           }
         */

         self->flow[ii] = -1.0*global_node_avg_density[ii]*global_node_avg_sep[ii]/global_node_rho_zero_density[ii];
         if (phi_ptr)
            self->flow[ii] -= global_node_avg_height[ii]*global_node_phi[ii]/global_node_rho_zero_density[ii];

         avgFlow += self->flow[ii];

      }


      //MPI_Allreduce(&avgFlow, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      tmp = avgFlow;
      avgFlow = tmp/(double)arraySize;

      free( node_avg_sep );
      free( node_avg_height );
      free( node_avg_density );
      free( node_rho_zero_density );
      if ( phi_ptr )
         free( node_phi );

      free( global_node_avg_sep );
      free( global_node_avg_height );
      free( global_node_avg_density );
      free( global_node_rho_zero_density );
      if ( phi_ptr )
         free( global_node_phi );

   }

   if (rank == 0)
   {
      printf("===========\n");
      printf("= Isostasy\n");
      printf("=\n");
      printf("= Average basal flow: %g\n", avgFlow);
      printf("===========\n");
   }
}


void lecode_tools_Isostasy_SetBC(int node_ind, int var_ind, void *_ctx, void* data, void *result)
{
   lecode_tools_Isostasy *self = lecode_tools_Isostasy_self;

   if ( self->avg )
   {
      ((double*)result)[0] = self->flow[0];
   }
   else
   {
      Grid *nodeGrid;
      int param[3], index, nDims;

      nodeGrid = *Mesh_GetExtension( self->mesh, Grid**,  self->mesh->vertGridId );
      Grid_Lift( nodeGrid, FeMesh_NodeDomainToGlobal( self->mesh, node_ind ), param );
      nDims = Mesh_GetDimSize( self->mesh );
      if ( nDims == 3 )
          index = param[0] + (param[self->zontalAxis]*nodeGrid->sizes[0]);
      else
          index = param[0];

      ((double*)result)[0] = self->flow[index];
   }
}


void lecode_tools_Isostasy_NonLinearEP(void* _sle, void* _ctx)
{
   lecode_tools_Isostasy *self = lecode_tools_Isostasy_self;
   SystemLinearEquations *sle = (SystemLinearEquations*)_sle;

   if ( self->maxIts == -1 || sle->nonLinearIteration_I < self->maxIts )
      lecode_tools_Isostasy_CalcBasalFlow( self );
}


void _lecode_tools_Isostasy_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data )
{
   lecode_tools_Isostasy* self = (lecode_tools_Isostasy*)component;
   ConditionFunction *cond_func;
   Dictionary* pluginDict = Codelet_GetPluginDictionary( component, cf->rootDict );

   lecode_tools_Isostasy_self = self;

   self->context = (AbstractContext*)Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"Context"  ), UnderworldContext, True, data );
   self->ctx = self->context;
   self->mesh = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"mesh"  ), FeMesh, True, data );
   self->swarm = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"swarm"  ), IntegrationPointsSwarm, True, data );
   /*
     self->heightField = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"HeightField"  ), FeVariable, True, data );
   */
   self->surfaceIdx = Stg_ComponentFactory_PluginGetInt( cf, self, "SurfaceIndex", -1 );
   self->vertAxis = Stg_ComponentFactory_PluginGetInt( cf, self, "VerticalAxis", 1 );
   if( self->vertAxis == 1 )
      self->zontalAxis = 2;
   else if( self->vertAxis == 2 )
      self->zontalAxis = 1;
   else
      Journal_Firewall( NULL,
       Journal_MyStream( Error_Type, self ),
       "\n\nError in %s for %s - Isostasy plugin requires 'VerticalAxis' to be '1' (standard UW) or '2'.",
       __func__,
       self->type );

   self->avg = Stg_ComponentFactory_PluginGetBool( cf, self, "UseAverage", True );
   self->vel_field = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"velocityField"  ), FeVariable, True, data );
   self->buoyancy = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"buoyancy"  ), BuoyancyForceTerm, False, data );
   self->rho_zero_mat = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"rhoZeroMaterial"  ), Material, True, data );
   self->sle = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"sle"  ), SystemLinearEquations, True, data );
   self->thermalSLE = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"thermalSLE"  ), SystemLinearEquations, False, data );
   if (self->thermalSLE)
   {
      self->tempField = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"temperatureField"  ), FeVariable, True, data );
      self->tempDotField = Stg_ComponentFactory_ConstructByName( cf, Dictionary_GetString( pluginDict, (Dictionary_Entry_Key)"temperatureDotField"  ), FeVariable, True, data );
   }
   self->maxIts = Stg_ComponentFactory_PluginGetInt( cf, self, "MaxIterations", -1 );

   SystemLinearEquations_SetToNonLinear(self->sle);
   SystemLinearEquations_AddNonLinearEP(self->sle, lecode_tools_Isostasy_Type, lecode_tools_Isostasy_NonLinearEP);

   cond_func = ConditionFunction_New( (ConditionFunction_ApplyFunc*)lecode_tools_Isostasy_SetBC, (Name)"lecode_tools_Isostasy", NULL );
   ConditionFunction_Register_Add( condFunc_Register, cond_func);

   /** check if PPCing */
   if(!self->buoyancy){
      self->ppcManager = NULL;
      self->ppcManager = Stg_ComponentFactory_PluginConstructByKey( cf, self, (Dictionary_Entry_Key)"Manager", PpcManager, False, data );
      if( !self->ppcManager  )
         self->ppcManager = Stg_ComponentFactory_ConstructByName( cf, (Name)"default_ppcManager", PpcManager, False, data  );
         if( !self->ppcManager )
            Journal_Firewall( NULL,
                 Journal_MyStream( Error_Type, self ),
                 "\n\nError in %s for %s - Neither a BuoyancyForceTerm or Ppc manager was found .\nIsostasy plugin requires one OR the other to operate.\n\n\n",
                 __func__,
                 self->type );
         if( self->thermalSLE )
            Journal_Firewall( NULL,
                 Journal_MyStream( Error_Type, self ),
                 "\n\nError in %s for %s - Ppc not compatible with thermalSLE option.\n\n\n",
                 __func__,
                 self->type );

      self->densityID = PpcManager_GetPpcFromDict( self->ppcManager, cf, self->name, (Dictionary_Entry_Key)"DensityLabel", "DensityLabel" );

   }

}


void lecode_tools_Isostasy_Build( void* component, void* data )
{
   lecode_tools_Isostasy *self = (lecode_tools_Isostasy*)component;
   Grid *nodeGrid;
   int arraySize, nDims;

   Stg_Component_Build( self->mesh, data, False );

   nDims = Mesh_GetDimSize( self->mesh );
   nodeGrid = *Mesh_GetExtension( self->mesh, Grid**,  self->mesh->vertGridId );

   arraySize=0;
   if ( nDims == 2 ) arraySize = nodeGrid->sizes[0];
   else if ( nDims == 3 ) arraySize = nodeGrid->sizes[0]*nodeGrid->sizes[self->zontalAxis];
   else assert(0);

   self->flow = (double*)malloc( arraySize*sizeof(double) );
   memset( self->flow, 0, arraySize*sizeof(double) );
}


void lecode_tools_Isostasy_Initialise( void* component, void* data ){}


void* _lecode_tools_Isostasy_DefaultNew( Name name )
{
   /* Variables set in this function */
   SizeT                                              _sizeOfSelf = sizeof(lecode_tools_Isostasy);
   Type                                                      type = lecode_tools_Isostasy_Type;
   Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
   Stg_Class_PrintFunction*                                _print = _Codelet_Print;
   Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
   Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lecode_tools_Isostasy_DefaultNew;
   Stg_Component_ConstructFunction*                    _construct = _lecode_tools_Isostasy_AssignFromXML;
   Stg_Component_BuildFunction*                            _build = lecode_tools_Isostasy_Build;
   Stg_Component_InitialiseFunction*                  _initialise = lecode_tools_Isostasy_Initialise;
   Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
   Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

   /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
   AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return _Codelet_New(  CODELET_PASSARGS  );
}


Index lecode_tools_Isostasy_Register( PluginsManager* pluginsManager )
{
   return PluginsManager_Submit( pluginsManager, lecode_tools_Isostasy_Type, (Name)"0", _lecode_tools_Isostasy_DefaultNew  );
}

