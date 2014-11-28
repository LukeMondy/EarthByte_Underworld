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

#ifndef __Q2ProjectionGenerator_h__
#define __Q2ProjectionGenerator_h__

	/** Textual name of this class */
	extern const Type Q2ProjectionGenerator_Type;

	/** Virtual function types */

	/** Q2ProjectionGenerator class contents */
	#define __Q2ProjectionGenerator			\
		/* General info */		\
		__ProjectionGenerator       \

	struct Q2ProjectionGenerator { __Q2ProjectionGenerator };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define C2SPHERICALGENERATOR_DEFARGS \
                CARTESIANGENERATOR_DEFARGS

	#define C2SPHERICALGENERATOR_PASSARGS \
                CARTESIANGENERATOR_PASSARGS

	void* _Q2ProjectionGenerator_DefaultNew( Name name );
	Q2ProjectionGenerator* _Q2ProjectionGenerator_New(  C2SPHERICALGENERATOR_DEFARGS  );
	void _Q2ProjectionGenerator_Init( Q2ProjectionGenerator* self );
	void _Q2ProjectionGenerator_Build( void* meshGenerator, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Q2ProjectionGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_Discretisaton_C2Generator_h__ */

