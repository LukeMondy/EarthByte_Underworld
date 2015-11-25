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

#ifndef __Spherical_Components_SphericalGenerator_h__
#define __Spherical_Components_SphericalGenerator_h__

/** Textual name of this class */
extern const Type SphericalGenerator_Type;

/** SphericalGenerator class contents */
#define __SphericalGenerator \
   /* General info */         \
   __CartesianGenerator       \
   /* note: the parametrisation - crdMin crdMax - */ \
   /* now represent spherical coordinate [theta, radius, phi] */ \
   double sph_res[3];  /* storage for the spherical resolution: order is [theta, radius, phi] */ \
   Bool fullAnnulus;  /* flag for FullAnnulus mesh */

struct SphericalGenerator
{
   __SphericalGenerator
};

/*--------------------------------------------------------------------------------------------------------------------------
** Constructors
*/




#ifndef ZERO
#define ZERO 0
#endif

#define SPHERICALGENERATOR_DEFARGS \
								CARTESIANGENERATOR_DEFARGS \
 
#define SPHERICALGENERATOR_PASSARGS \
								CARTESIANGENERATOR_PASSARGS \
 
void* _SphericalGenerator_DefaultNew( Name name );
SphericalGenerator* _SphericalGenerator_New(  SPHERICALGENERATOR_DEFARGS  );

void _SphericalGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
void _SphericalGenerator_Initialise( void* meshGenerator, void* data );

void _SphericalGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh );
void SphericalGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );


#endif /* __StgDomain_Mesh_SphericalGenerator_h__ */

