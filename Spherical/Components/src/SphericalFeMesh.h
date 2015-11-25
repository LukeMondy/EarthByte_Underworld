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

#ifndef __SphericalFeMesh_h__
#define __SphericalFeMesh_h__

	/** Textual name of this class */
	extern const Type SphericalFeMesh_Type;

	/** Virtual function types */

	/** Class contents */
	#define __SphericalFeMesh				\
		/* General info */			\
		__FeMesh					\
							\
    IndexSet*           innerSet;  /* inner radius bnd nodes */ \
    IndexSet*           outerSet; /* outer radius bnd nodes */  \
    IndexSet*           innerElSet;  /* inner radius bnd Els */ \
    IndexSet*           outerElSet; /* outer radius bnd Els */  \
		

	struct SphericalFeMesh { __SphericalFeMesh };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SPHERICALFEMESH_DEFARGS \
                FEMESH_DEFARGS

	#define SPHERICALFEMESH_PASSARGS \
                FEMESH_PASSARGS

	SphericalFeMesh* SphericalFeMesh_New( Name name, AbstractContext* context );
	SphericalFeMesh* _SphericalFeMesh_New(  FEMESH_DEFARGS  );
	void _SphericalFeMesh_Init( SphericalFeMesh* self, ElementType* elType, const char* family, Bool elementMesh );
   SphericalFeMesh* _SphericalFeMesh_DefaultNew( Name name );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _SphericalFeMesh_Delete( void* feMesh );
	void _SphericalFeMesh_Print( void* feMesh, Stream* stream );
	void _SphericalFeMesh_AssignFromXML( void* feMesh, Stg_ComponentFactory* cf, void* data );
	void _SphericalFeMesh_Build( void* feMesh, void* data );
	void _SphericalFeMesh_Initialise( void* feMesh, void* data );
	void _SphericalFeMesh_Execute( void* feMesh, void* data );
	void _SphericalFeMesh_Destroy( void* feMesh, void* data );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void SphericalFeMesh_SetElementFamily( void* feMesh, const char* family );
	void SphericalFeMesh_SetElementType( void* feMesh, ElementType* elType );

	ElementType* SphericalFeMesh_GetElementType( void* feMesh, unsigned element );

	unsigned SphericalFeMesh_GetNodeLocalSize( void* feMesh );
	unsigned SphericalFeMesh_GetNodeRemoteSize( void* feMesh );
	unsigned SphericalFeMesh_GetNodeDomainSize( void* feMesh );
	unsigned SphericalFeMesh_GetNodeGlobalSize( void* feMesh );
	unsigned SphericalFeMesh_GetElementLocalSize( void* feMesh );
	unsigned SphericalFeMesh_GetElementDomainSize( void* feMesh );
	unsigned SphericalFeMesh_GetElementRemoteSize( void* feMesh );
	unsigned SphericalFeMesh_GetElementGlobalSize( void* feMesh );

	unsigned SphericalFeMesh_GetElementNodeSize( void* feMesh, unsigned element );
	unsigned SphericalFeMesh_GetNodeElementSize( void* feMesh, unsigned node );
	void SphericalFeMesh_GetElementNodes( void* feMesh, unsigned element, IArray* inc );
	void SphericalFeMesh_GetNodeElements( void* feMesh, unsigned node, IArray* inc );

	unsigned SphericalFeMesh_ElementDomainToGlobal( void* feMesh, unsigned domain );
	Bool SphericalFeMesh_ElementGlobalToDomain( void* feMesh, unsigned global, unsigned* domain );
	unsigned SphericalFeMesh_NodeDomainToGlobal( void* feMesh, unsigned domain );
	Bool SphericalFeMesh_NodeGlobalToDomain( void* feMesh, unsigned global, unsigned* domain );

	void SphericalFeMesh_CoordGlobalToLocal( void* feMesh, unsigned element, double* global, double* local );
	void SphericalFeMesh_CoordLocalToGlobal( void* feMesh, unsigned element, double* local, double* global );
	void SphericalFeMesh_EvalBasis( void* feMesh, unsigned element, double* localCoord, double* basis );
	void SphericalFeMesh_EvalLocalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs );
	void SphericalFeMesh_EvalGlobalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs, double* jacDet );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void SphericalFeMesh_Destruct( SphericalFeMesh* self );

#endif /* __StgFEM_Discretisaton_SphericalFeMesh_h__ */

