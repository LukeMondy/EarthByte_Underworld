/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	AuScope - http://www.auscope.org
**	Monash University, AuScope SAM VIC - http://www.auscope.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
**	Abstract class defining the interface for a System of Linear Equations solver.
**
** Assumptions:
**
** Comments:
**	Note - as 1 September 2004 (rev 1994), the functioality for building Matrices etc
**	that was in this class is back in the SLIntegrator_Polar class.
**
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Spherical_Components_SLIntegrator_Polar_h__
#define __Spherical_Components_SLIntegrator_Polar_h__

   /** Textual name of this class */
   extern const Type SLIntegrator_Polar_Type;

   /** SLIntegrator_Polar class contents */
   #define __SLIntegrator_Polar \
      __Stg_Component \
      \
      /* SLIntegrator_Polar info */ \
      FeVariable*                  velocityField;	   \
      Stg_ObjectList*              variableList;       	   \
      Stg_ObjectList*              varStarList;  	   \
      Bool*                        pureAdvection;          \
      FiniteElementContext*        context;                \
      double*			   abcissa;		   \
      double*			   Ni;			   \
      double**			   GNix;		   \
      IArray*			   inc;			   \
      double			   courant;		   \
      unsigned**		   elPatch;		   \

   /** Abstract class defining the interface for a SLIntegrator_Polar solver - see SLIntegrator_Polar.h */
   struct SLIntegrator_Polar { __SLIntegrator_Polar };

   /* --- Constructor functions --- */
   void* _SLIntegrator_Polar_DefaultNew( Name name );

   /** Creation implementation */

   #ifndef ZERO
   #define ZERO 0
   #endif

   #define SLINTEGRATOR_POLAR_DEFARGS \
                STG_COMPONENT_DEFARGS

   #define SLINTEGRATOR_POLAR_PASSARGS \
                STG_COMPONENT_PASSARGS

   SLIntegrator_Polar* _SLIntegrator_Polar_New(  SLINTEGRATOR_POLAR_DEFARGS  );

   /* --- Virtual function implementations --- */

   /** Class Virtual Functions Implementations */
   void* _SLIntegrator_Polar_Copy( void* slIntegrator, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
   void _SLIntegrator_Polar_Delete( void* slIntegrator );
   void _SLIntegrator_Polar_Print( void* slIntegrator, Stream* stream );

   /** Stg_Component_Build() implementation: does nothing by default as some solvers may not need it. */
   void _SLIntegrator_Polar_Build( void* slIntegrator, void* data );

   void _SLIntegrator_Polar_AssignFromXML( void* slIntegrator, Stg_ComponentFactory* cf, void* data );

   /** Stg_Component_Initialise() implementation: does nothing by default as some solvers may not neet it. */
   void _SLIntegrator_Polar_Initialise( void* slIntegrator, void* data );

   /** Stg_Component_Execute() implementation: Calls SolverSetup() for any per-iteration setup, then
   calls Solve() to calculate the new solutions. */
   void _SLIntegrator_Polar_Execute( void* slIntegrator, void* data );

   void _SLIntegrator_Polar_Destroy( void* slIntegrator, void* data );

   /* --- Private Functions --- */
   void SLIntegrator_Polar_InitSolve( void* slIntegrator, void* data );

   /* --- Public functions --- */
   void SLIntegrator_Polar_CubicInterpolator( void* slIntegrator, FeVariable* feVariable, double* coord, double* result );
   void SLIntegrator_Polar_PeriodicUpdate( double* pos, double* min, double* max, Bool* isPeriodic );
   void SLIntegrator_Polar_CheckBounds( double* pos, double min, double max );
   void SLIntegrator_Polar_IntegrateRK4( void* slIntegrator, FeVariable* velocityField, double dt, double* origin, double* position );

   /** Solve:- calculate the new values for all solution vectors in the system. */
   void SLIntegrator_Polar_Solve( void* slIntegrator, FeVariable* variableField, FeVariable* variableFieldPrime );

   void SLIntegrator_Polar_ShapeFuncs( void* slIntegrator, const double lCoord[], double* const Ni );
   void SLIntegrator_Polar_ShapeFuncDerivs( void* slIntegrator, const double lCoord[], double** const GNix );
   void SLIntegrator_Polar_GlobalToLocal( void* slIntegrator, void* _mesh, unsigned* nodeInds, const double* gCoord, double* lCoord );

   double SLIntegrator_Polar_CalcAdvDiffDt( void* slIntegrator, FiniteElementContext* context );

   void SLIntegrator_Polar_InitPatches( void* slIntegrator );

#endif

