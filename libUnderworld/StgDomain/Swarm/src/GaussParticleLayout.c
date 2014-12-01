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
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "shortcuts.h"
#include "ParticleLayout.h"
#include "PerCellParticleLayout.h"
#include "GaussParticleLayout.h"

#include "SwarmClass.h"
#include "StandardParticle.h"
#include "IntegrationPoint.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

const Type GaussParticleLayout_Type = "GaussParticleLayout";

GaussParticleLayout* GaussParticleLayout_New( 
   Name name, 
   AbstractContext* context,
   CoordSystem      coordSystem,
   Bool             weightsInitialisedAtStartup,
   Dimension_Index dim, 
   Particle_InCellIndex* particlesPerDim ) {

   GaussParticleLayout* self = _GaussParticleLayout_DefaultNew( name );
   _ParticleLayout_Init( self, context, coordSystem, weightsInitialisedAtStartup );
   _PerCellParticleLayout_Init( self );
   _GaussParticleLayout_Init( self, dim, particlesPerDim );

   return self;
}

GaussParticleLayout* _GaussParticleLayout_New(  GAUSSPARTICLELAYOUT_DEFARGS  )
{
	GaussParticleLayout* self;
	
   /* hard-wire these */
   coordSystem = LocalCoordSystem;
   weightsInitialisedAtStartup = True;
   nameAllocationType = NON_GLOBAL;

	/* Allocate memory */
	self = (GaussParticleLayout*)_PerCellParticleLayout_New(  PERCELLPARTICLELAYOUT_PASSARGS  );

   self->dim = dim;
   if( particlesPerDim )
     memcpy( self->particlesPerDim, particlesPerDim, 3 * sizeof(Particle_InCellIndex) );
	
	return self;
}

void _GaussParticleLayout_Init( void* gaussParticleLayout, Dimension_Index dim, Particle_InCellIndex* particlesPerDim ) {
   GaussParticleLayout* self = (GaussParticleLayout*)gaussParticleLayout;

   self->isConstructed       = True;
   self->coordSystem         = LocalCoordSystem;
   self->weightsInitialisedAtStartup = True;
   self->dim                 = dim;
   memcpy( self->particlesPerDim, particlesPerDim, 3 * sizeof(Particle_InCellIndex) );
}

void _GaussParticleLayout_Delete( void* gaussParticleLayout ) {
	GaussParticleLayout* self = (GaussParticleLayout*)gaussParticleLayout;
	
	_PerCellParticleLayout_Delete( self );
}

void _GaussParticleLayout_Print( void* gaussParticleLayout, Stream* stream ) {
	GaussParticleLayout* self = (GaussParticleLayout*)gaussParticleLayout;
	
	/* General info */
	Journal_Printf( stream, "GaussParticleLayout (ptr): %p:\n", self );
	
	/* Parent class info */
	_PerCellParticleLayout_Print( self, stream );
	
	/* Virtual info */
	
	/* GaussParticleLayout */
	Stream_Indent( stream );
	Journal_PrintValue( stream, self->dim );
	Journal_PrintArray( stream, self->particlesPerDim, self->dim );
	Stream_UnIndent( stream );
}


void* _GaussParticleLayout_Copy( void* gaussParticleLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	GaussParticleLayout*	self = (GaussParticleLayout*)gaussParticleLayout;
	GaussParticleLayout*	newGaussParticleLayout;
	
	newGaussParticleLayout = (GaussParticleLayout*)_PerCellParticleLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	newGaussParticleLayout->dim = self->dim;
	memcpy( newGaussParticleLayout->particlesPerDim, self->particlesPerDim, 3 * sizeof(unsigned int) );

	return (void*)newGaussParticleLayout;
}

void* _GaussParticleLayout_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                                     _sizeOfSelf = sizeof(GaussParticleLayout);
	Type                                                                             type = GaussParticleLayout_Type;
	Stg_Class_DeleteFunction*                                                     _delete = _GaussParticleLayout_Delete;
	Stg_Class_PrintFunction*                                                       _print = _GaussParticleLayout_Print;
	Stg_Class_CopyFunction*                                                         _copy = _GaussParticleLayout_Copy;
	Stg_Component_DefaultConstructorFunction*                         _defaultConstructor = _GaussParticleLayout_DefaultNew;
	Stg_Component_ConstructFunction*                                           _construct = _GaussParticleLayout_AssignFromXML;
	Stg_Component_BuildFunction*                                                   _build = _GaussParticleLayout_Build;
	Stg_Component_InitialiseFunction*                                         _initialise = _GaussParticleLayout_Initialise;
	Stg_Component_ExecuteFunction*                                               _execute = _GaussParticleLayout_Execute;
	Stg_Component_DestroyFunction*                                               _destroy = _GaussParticleLayout_Destroy;
	AllocationType                                                     nameAllocationType = NON_GLOBAL;
	ParticleLayout_SetInitialCountsFunction*                            _setInitialCounts = _PerCellParticleLayout_SetInitialCounts;
	ParticleLayout_InitialiseParticlesFunction*                      _initialiseParticles = _PerCellParticleLayout_InitialiseParticles;
	CoordSystem                                                               coordSystem = LocalCoordSystem;
	Bool                                                      weightsInitialisedAtStartup = True;
	PerCellParticleLayout_InitialCountFunction*                             _initialCount = _GaussParticleLayout_InitialCount;
	PerCellParticleLayout_InitialiseParticlesOfCellFunction*   _initialiseParticlesOfCell = _GaussParticleLayout_InitialiseParticlesOfCell;
	Dimension_Index                                                                   dim = 0;
	Particle_InCellIndex*                                                 particlesPerDim = NULL;

	return (GaussParticleLayout*)_GaussParticleLayout_New(  GAUSSPARTICLELAYOUT_PASSARGS  );
}

void _GaussParticleLayout_AssignFromXML( void* gaussParticleLayout, Stg_ComponentFactory* cf, void* data ) {
	GaussParticleLayout*   self = (GaussParticleLayout*)gaussParticleLayout;
	Particle_InCellIndex   particlesPerDim[3];
	Particle_InCellIndex   defaultVal;
	Dimension_Index        dim;

   _PerCellParticleLayout_AssignFromXML( self, cf, data );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, (Dictionary_Entry_Key)"dim", 0  );

	defaultVal = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"gaussParticles", 2  );

	particlesPerDim[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"gaussParticlesX", defaultVal  );
	particlesPerDim[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"gaussParticlesY", defaultVal );
	if ( dim == 3  )
		particlesPerDim[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"gaussParticlesZ", defaultVal  );
	else
		particlesPerDim[ K_AXIS ] = 1;	

	_GaussParticleLayout_Init( self, dim, particlesPerDim );
}
	
void _GaussParticleLayout_Build( void* gaussParticleLayout, void* data ) {
}
	
void _GaussParticleLayout_Initialise( void* gaussParticleLayout, void* data ) {
}
	
void _GaussParticleLayout_Execute( void* gaussParticleLayout, void* data ) {
}

void _GaussParticleLayout_Destroy( void* gaussParticleLayout, void* data ) {
}

Particle_InCellIndex _GaussParticleLayout_InitialCount( void* gaussParticleLayout, void* celllayout, Cell_Index cell_I )
{
	GaussParticleLayout* self   = (GaussParticleLayout*)gaussParticleLayout;
	Particle_InCellIndex count;
	Dimension_Index      dim;	
	Dimension_Index      dim_I;
	
	dim = self->dim;
	count = 1;
	for( dim_I = 0; dim_I < dim; dim_I++ ) {
		count = count * (Particle_InCellIndex)( self->particlesPerDim[ dim_I ] );
	}
	
	return count;
	
}

#define TRIPLE_MAX( A, B, C )  MAX( MAX( (A), (B) ), (C) )

/* remember this only has to initialise one particle at a time */
void _GaussParticleLayout_InitialiseParticlesOfCell( void* gaussParticleLayout, void* _swarm, Cell_Index cell_I )
{
	GaussParticleLayout*      self                = (GaussParticleLayout*)gaussParticleLayout;
	Swarm*                    swarm               = (Swarm*)_swarm;
	IntegrationPoint*         particle            = NULL;
	Index                     index2D;
	Particle_InCellIndex      maxParticlesPerDim;
	IJK                       ijkIndex;
	Index                     index;
	Dimension_Index           dim_I;
	div_t                     divide;
	double*                   weights;
	double*                   abscissa;
	Coord                     min;
	Coord                     max;
	Particle_InCellIndex      particlesThisCell = swarm->cellParticleCountTbl[cell_I];
	Particle_InCellIndex      cParticle_I = 0;
	

	if ( 0 == strcmp( swarm->type, "MaterialPointsSwarm" ) ) {
		/* TODO: This is a special rule allowing a Gauss particle layout to be used to initialise
		global co-ords if you want to use it in a material swarm */
		self->coordSystem = GlobalCoordSystem;
	}

	Swarm_GetCellMinMaxCoords( _swarm, cell_I, min, max );

	/* Allocate Memory */
	maxParticlesPerDim = TRIPLE_MAX( self->particlesPerDim[ I_AXIS ],
		self->particlesPerDim[ J_AXIS ], self->particlesPerDim[ K_AXIS ] );

	weights   = Memory_Alloc_Array( double, maxParticlesPerDim, "gauss weights" );
	abscissa  = Memory_Alloc_Array( double, maxParticlesPerDim, "gauss abscissa" );

	for ( cParticle_I = 0; cParticle_I < particlesThisCell; cParticle_I++ ) {
		particle = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		particle->owningCell = cell_I;
		
		/* Find the i, j, k index of this particular particle */
		divide = div( cParticle_I, self->particlesPerDim[ I_AXIS ] * self->particlesPerDim[ J_AXIS ] );
		ijkIndex[ K_AXIS ] = divide.quot;
		index2D = divide.rem;

		divide = div( index2D, self->particlesPerDim[ I_AXIS ] );
		ijkIndex[ J_AXIS ] = divide.quot;
		ijkIndex[ I_AXIS ] = divide.rem;

		particle->weight = 1.0;
		for( dim_I = 0 ; dim_I < self->dim ; dim_I++ ) {
			index = ijkIndex[ dim_I ];
			GaussParticleLayout_GetAbscissaAndWeights1D( weights, abscissa, self->particlesPerDim[ dim_I ] );

			/* Assign particle stuff */
			/* TODO: this is a hack, since this class doesn't officially know that the MaterialPointsSwarm
			 * exists yet. However, we need some way for material swarms to use this layout, for testing
			 * purposes. In the simple system of only 1 swarm type, this component always initialised
			 * both local and global co-ords.
			 * -- PatrickSunter - 1 May 2006
			 */
			if ( 0 == strcmp( swarm->type, "MaterialPointsSwarm" ) ) {
				((GlobalParticle*)particle)->coord[ dim_I ] = 
					min[ dim_I ] +
						0.5 * ( max[ dim_I ] - min[ dim_I ] ) 
						* ( abscissa[ index ] + 1.0 );
			}
			else {
				particle->xi[ dim_I ] = abscissa[ index ];
				particle->weight *= weights[ index ];
			}	
		}
		
	}	

	Memory_Free( weights );
	Memory_Free( abscissa );
}


/* Values taken from table from:
 * Eric W. Weisstein. "Legendre-Gauss Quadrature." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Legendre-GaussQuadrature.html */
/* (Mirko): Values for >5 calculated by following Maple code:
restart:
kernelopts(printbytes=false);
kernelopts(datalimit=200000);
f:=(n,x)->simplify(LegendreP(n,x),'LegendreP'):
Digits:=45:
fname:=cat("gauss_weight",Digits):
for n from 6 to 11 do
 sols:=fsolve(f(n,x)=0,x);
 fd := fopen(fname,APPEND):
 fprintf(fd,"case %d:\n",n):
 print(n);
 for i from 1 to n do
   a:=Re( evalf(sols[i],Digits) ):w:=Re( evalf(2*(1-a^2)/( ((n+1)^2)*(f(n+1,a))^2 ),Digits) ):
   print(a,w);
   fprintf(fd,"abscissa[%d] = %.25f; weight[%d] =%.25f;\n",i-1,a,i-1,w):
 od:
 fprintf(fd,"break;\n"):
 fclose(fd):
od:

*/
void GaussParticleLayout_GetAbscissaAndWeights1D( double* weight, double* abscissa, Index pointCount ) {
  switch ( pointCount ) {
  case 1:
    abscissa[0]  = 0.0;
    weight[0] = 2.0;
    break;
  case 2:
    abscissa[0]  = - 1.0/sqrt(3.0);
    abscissa[1]  = - abscissa[0];
		
    weight[0] = 1.0;
    weight[1] = weight[0];
    break;
  case 3:
    abscissa[0]  = - sqrt(15.0)/5.0;
    abscissa[1]  = 0.0;
    abscissa[2]  = - abscissa[0];
		
    weight[0] = 5.0/9.0;
    weight[1] = 8.0/9.0;
    weight[2] = weight[0];
    break;
  case 4:
    abscissa[0]  = - sqrt( 525.0 + 70.0 * sqrt(30.0) )/35.0;
    abscissa[1]  = - sqrt( 525.0 - 70.0 * sqrt(30.0) )/35.0;
    abscissa[2]  = - abscissa[1];
    abscissa[3]  = - abscissa[0];
			
    weight[0] = (18.0 - sqrt(30.0))/36.0;
    weight[1] = (18.0 + sqrt(30.0))/36.0;
    weight[2] = weight[1];
    weight[3] = weight[0];
    break;
  case 5:
    abscissa[0]  = - sqrt( 245.0 + 14.0 * sqrt( 70.0 ) )/21.0;
    abscissa[1]  = - sqrt( 245.0 - 14.0 * sqrt( 70.0 ) )/21.0;
    abscissa[2]  = 0.0;
    abscissa[3]  = - abscissa[1];
    abscissa[4]  = - abscissa[0];

    weight[0] = ( 322.0 - 13.0 * sqrt( 70.0 ) )/900.0;
    weight[1] = ( 322.0 + 13.0 * sqrt( 70.0 ) )/900.0;
    weight[2] = 128.0/225.0;
    weight[3] = weight[1];
    weight[4] = weight[0];
    break;
  case 6:
    abscissa[0] = -0.9324695142031520278123016; weight[0] =0.1713244923791703450402961;
    abscissa[1] = -0.6612093864662645136613996; weight[1] =0.3607615730481386075698335;
    abscissa[2] = -0.2386191860831969086305017; weight[2] =0.4679139345726910473898703;
    abscissa[3] = 0.2386191860831969086305017;  weight[3] =0.4679139345726910473898703;
    abscissa[4] = 0.6612093864662645136613996;  weight[4] =0.3607615730481386075698335;
    abscissa[5] = 0.9324695142031520278123016;  weight[5] =0.1713244923791703450402961;
    break;
  case 7:
    abscissa[0] = -0.9491079123427585245261897; weight[0] =0.1294849661688696932706114;
    abscissa[1] = -0.7415311855993944398638648; weight[1] =0.2797053914892766679014678;
    abscissa[2] = -0.4058451513773971669066064; weight[2] =0.3818300505051189449503698;
    abscissa[3] =  0.0000000000000000000000000; weight[3] =0.4179591836734693877551020;
    abscissa[4] =  0.4058451513773971669066064; weight[4] =0.3818300505051189449503698;
    abscissa[5] =  0.7415311855993944398638648; weight[5] =0.2797053914892766679014678;
    abscissa[6] =  0.9491079123427585245261897; weight[6] =0.1294849661688696932706114;
    break;
  case 8:
    abscissa[0] = -0.9602898564975362316835609; weight[0] =0.1012285362903762591525314;
    abscissa[1] = -0.7966664774136267395915539; weight[1] =0.2223810344533744705443560;
    abscissa[2] = -0.5255324099163289858177390; weight[2] =0.3137066458778872873379622;
    abscissa[3] = -0.1834346424956498049394761; weight[3] =0.3626837833783619829651504;
    abscissa[4] =  0.1834346424956498049394761; weight[4] =0.3626837833783619829651504;
    abscissa[5] =  0.5255324099163289858177390; weight[5] =0.3137066458778872873379622;
    abscissa[6] =  0.7966664774136267395915539; weight[6] =0.2223810344533744705443560;
    abscissa[7] =  0.9602898564975362316835609; weight[7] =0.1012285362903762591525314;
    break;
  case 9:
    abscissa[0] = -0.9681602395076260898355762; weight[0] =0.0812743883615744119718922;
    abscissa[1] = -0.8360311073266357942994298; weight[1] =0.1806481606948574040584720;
    abscissa[2] = -0.6133714327005903973087020; weight[2] =0.2606106964029354623187429;
    abscissa[3] = -0.3242534234038089290385380; weight[3] =0.3123470770400028400686304;
    abscissa[4] = 0.0000000000000000000000000; weight[4] =0.3302393550012597631645251;
    abscissa[5] = 0.3242534234038089290385380; weight[5] =0.3123470770400028400686304;
    abscissa[6] = 0.6133714327005903973087020; weight[6] =0.2606106964029354623187429;
    abscissa[7] = 0.8360311073266357942994298; weight[7] =0.1806481606948574040584720;
    abscissa[8] = 0.9681602395076260898355762; weight[8] =0.0812743883615744119718922;
    break;
  case 10:
    abscissa[0] = -0.9739065285171717200779640; weight[0] =0.0666713443086881375935688;
    abscissa[1] = -0.8650633666889845107320967; weight[1] =0.1494513491505805931457763;
    abscissa[2] = -0.6794095682990244062343274; weight[2] =0.2190863625159820439955349;
    abscissa[3] = -0.4333953941292471907992659; weight[3] =0.2692667193099963550912269;
    abscissa[4] = -0.1488743389816312108848260; weight[4] =0.2955242247147528701738930;
    abscissa[5] = 0.1488743389816312108848260; weight[5] =0.2955242247147528701738930;
    abscissa[6] = 0.4333953941292471907992659; weight[6] =0.2692667193099963550912269;
    abscissa[7] = 0.6794095682990244062343274; weight[7] =0.2190863625159820439955349;
    abscissa[8] = 0.8650633666889845107320967; weight[8] =0.1494513491505805931457763;
    abscissa[9] = 0.9739065285171717200779640; weight[9] =0.0666713443086881375935688;
    break;
  case 11:
    abscissa[0] = -0.9782286581460569928039380; weight[0] =0.0556685671161736664827537;
    abscissa[1] = -0.8870625997680952990751578; weight[1] =0.1255803694649046246346943;
    abscissa[2] = -0.7301520055740493240934163; weight[2] =0.1862902109277342514260976;
    abscissa[3] = -0.5190961292068118159257257; weight[3] =0.2331937645919904799185237;
    abscissa[4] = -0.2695431559523449723315320; weight[4] =0.2628045445102466621806889;
    abscissa[5] = 0.0000000000000000000000000; weight[5] =0.2729250867779006307144835;
    abscissa[6] = 0.2695431559523449723315320; weight[6] =0.2628045445102466621806889;
    abscissa[7] = 0.5190961292068118159257257; weight[7] =0.2331937645919904799185237;
    abscissa[8] = 0.7301520055740493240934163; weight[8] =0.1862902109277342514260976;
    abscissa[9] = 0.8870625997680952990751578; weight[9] =0.1255803694649046246346943;
    abscissa[10] = 0.9782286581460569928039380; weight[10] =0.0556685671161736664827537;
    break;
  default:
    Journal_Firewall( 
      pointCount <= 11, 
      Journal_Register( Error_Type, (Name)GaussParticleLayout_Type  ),
      "In func %s: Cannot give values for '%u' gauss points\n.", 
      __func__, 
      pointCount );
    exit(1);
  }
}


