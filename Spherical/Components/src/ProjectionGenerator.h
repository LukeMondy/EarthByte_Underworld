/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Spherical_Components_ProjectionGenerator_h__
#define __Spherical_Components_ProjectionGenerator_h__

/** Textual name of this class */
extern const Type ProjectionGenerator_Type;

typedef double (RingMappingFunc) ( void* _self, int sq_i );
/** ProjectionGenerator class contents */
#define __ProjectionGenerator \
   __CartesianGenerator       \
   RingMappingFunc* _mapRing;  \
   int *nodeSqId; /* stores which cube/square the node belongs to in projection algorithm */ \
   int nSquares;              \
   double distance;     \
   double superel; \
   double decay;        \
   Bool equiangle; \
   Bool rotateMesh; \
   double pgen_res[3];


void _ProjectionGenerator_Destroy( void* meshGenerator, void* data );

struct ProjectionGenerator
{
  __ProjectionGenerator
};

/*--------------------------------------------------------------------------------------------------------------------------
** Constructors
*/




#ifndef ZERO
#define ZERO 0
#endif

#define PROJECTIONGENERATOR_DEFARGS \
								CARTESIANGENERATOR_DEFARGS \
 
#define PROJECTIONGENERATOR_PASSARGS \
								CARTESIANGENERATOR_PASSARGS \
 
ProjectionGenerator* _ProjectionGenerator_DefaultNew( Name name );
ProjectionGenerator* _ProjectionGenerator_New(  PROJECTIONGENERATOR_DEFARGS  );

void _ProjectionGenerator_Init( ProjectionGenerator* self, int nSquares, int map, double distance, double decay, double superel );

void ProjectionGenerator_Generate( void* meshGenerator, void* data );
void _ProjectionGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
void _ProjectionGenerator_Build( void* meshGenerator, void* data );
void _ProjectionGenerator_Initialise( void* meshGenerator, void* data );

void _ProjectionGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh );
void ProjectionGenerator_EquiDistanceGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );

void ProjectionGenerator_EquiAngleGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );
void ProjectionGenerator_EquiAngle3DGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );

void ProjectionGenerator_CalcFullAnnulus( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );

double ProjectionGenerator_PowerMap( void* _self, int sq_i );
double ProjectionGenerator_ExpMap( void* _self, int sq_i );



#endif /* __StgDomain_Mesh_ProjectionGenerator_h__ */

