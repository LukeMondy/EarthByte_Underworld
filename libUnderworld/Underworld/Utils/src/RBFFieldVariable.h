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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Utils_RBFFieldVariable_h__
#define __Underworld_Utils_RBFFieldVariable_h__

	/** Textual name of this class */
	extern const Type RBFFieldVariable_Type;
	
	/** RBFFieldVariable contents */
	#define __RBFFieldVariable                                 \
		/* General info */                                      \
		__FieldVariable                                         \
		                                                        \
		/* Member info */                                       \
      RBFManager*                 rbfManager;                 \
		Bool                        useShepardCorrection;       \
		double                      defaultValue;               \
		double                      offset;                     \
		double                      shepardCorrectionThreshold; \
      SwarmVariable*              swarmVariable;

	struct RBFFieldVariable { __RBFFieldVariable };	



	/** Creation implementation */
	void* _RBFFieldVariable_DefaultNew( Name name );

	
	#ifndef ZERO
	#define ZERO 0
	#endif

   #define RBFFIELDVARIABLE_DEFARGS \
                FIELDVARIABLE_DEFARGS
                
   #define RBFFIELDVARIABLE_PASSARGS \
                FIELDVARIABLE_PASSARGS  

   RBFFieldVariable* _RBFFieldVariable_New(  RBFFIELDVARIABLE_DEFARGS  );

	/** Member initialisation implementation */
	void _RBFFieldVariable_Init( 
   RBFFieldVariable*          self,
   RBFManager*                rbfManager,
   Bool                       useShepardCorrection,
   double                     defaultValue,
   double                     offset,
   double                     shepardCorrectionThreshold,
   SwarmVariable*             swarmVariable ) ;
	
	void _RBFFieldVariable_AssignFromXML( void* RBFFieldVariable, Stg_ComponentFactory* cf, void* data ) ;

	void _RBFFieldVariable_Build( void* RBFFieldVariable, void* data ) ;

	void _RBFFieldVariable_Execute( void* RBFFieldVariable, void* data ) ;

	void _RBFFieldVariable_Destroy( void* RBFFieldVariable, void* data ) ;

	void _RBFFieldVariable_Initialise( void* RBFFieldVariable, void* data ) ;
	
   /** Interpolate the value of the RBF variable at a particular coord **/
	InterpolationResult _RBFFieldVariable_InterpolateValueAt( void* RBFFieldVariable, double* coord, double* value );

	/** Implementations of the min and max val functions */
	void _RBFFieldVariable_CacheMinMaxGlobalFieldMagnitude( void* feVariable );

	InterpolationResult  RBFFieldVariable_InterpolateGradientValueAt( void* _RBFFieldVariable, double* coord, double* value, int axis );



#endif /* __Underworld_Utils_RBFFieldVariable_h__ */

