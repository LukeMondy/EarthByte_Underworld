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

#ifndef __Spherical_Components_Q2RSGenerator_h__
#define __Spherical_Components_Q2RSGenerator_h__

/** Textual name of this class */
extern const Type Q2RSGenerator_Type;

/** Q2RSGenerator class contents */
#define __Q2RSGenerator \
   /* General info */         \
   __RSGenerator       \

struct Q2RSGenerator
{
   __Q2RSGenerator
};

/*--------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

#ifndef ZERO
#define ZERO 0
#endif

#define Q2RSGENERATOR_DEFARGS \
								RSGENERATOR_DEFARGS \
 
#define Q2RSGENERATOR_PASSARGS \
								RSGENERATOR_PASSARGS \
 
void* _Q2RSGenerator_DefaultNew( Name name, AbstractContext* context );
Q2RSGenerator* _Q2RSGenerator_New( Q2RSGENERATOR_DEFARGS );

void _Q2RSGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );
void _Q2RSGenerator_Initialise( void* meshGenerator, void* data );

void _Q2RSGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh );
void Q2RSGenerator_CalcCurvilinearGeom( void* _self, Mesh* mesh, Sync* sync, Grid* grid, unsigned* inds, double* steps );


#endif /* __StgDomain_Mesh_Q2RSGenerator_h__ */

