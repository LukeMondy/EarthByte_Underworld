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

#ifndef __Spherical_Components_RSGenerator_h__
#define __Spherical_Components_RSGenerator_h__

/** Textual name of this class */
extern const Type RSGenerator_Type;

/** RSGenerator class contents */
#define __RSGenerator \
   /* General info */         \
   __CartesianGenerator       \
   /* note: the parametrisation - crdMin crdMax - */ \
   /* now represent spherical coordinate [theta, radius, phi] */ \
   double sph_res[3];  /* storage for the spherical resolution: order is [theta, radius, phi] */ \
   Bool fullAnnulus;  /* flag for FullAnnulus mesh */ \
   Bool sixthOnly;  /* flag for 3D sixth of a sphere */

struct RSGenerator
{
   __RSGenerator
};

/*--------------------------------------------------------------------------------------------------------------------------
** Constructors
*/




#ifndef ZERO
#define ZERO 0
#endif

#define RSGENERATOR_DEFARGS \
								CARTESIANGENERATOR_DEFARGS \
 
#define RSGENERATOR_PASSARGS \
								CARTESIANGENERATOR_PASSARGS \
 
void* _RSGenerator_DefaultNew( Name name );
RSGenerator* _RSGenerator_New(  RSGENERATOR_DEFARGS  );

void _RSGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
void _RSGenerator_Initialise( void* meshGenerator, void* data );

void _RSGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh );
void RSGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );


#endif /* __StgDomain_Mesh_RSGenerator_h__ */

