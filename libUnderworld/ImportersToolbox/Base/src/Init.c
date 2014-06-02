/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005-2010, Monash University
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice,
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright
**			notice, this list of conditions and the following disclaimer in the
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors
**			may be used to endorse or promote products derived from this software
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%  Louis Moresi - Louis.Moresi(at)monash.edu
*%
*% Development Team :
*%  http://www.underworldproject.org/aboutus.html
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "Base.h"

Stream* ImportersToolboxInfo  = NULL;
Stream* ImportersToolboxDebug = NULL;
Stream* ImportersToolboxError = NULL;

Bool ImportersToolbox_Base_Init()
{
   Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

   Journal_Printf( Journal_Register( Debug_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

   /* Setup some streams for you to use */
   ImportersToolboxInfo  = Journal_Register(  Info_Type, (Name)"ImportersToolboxInfo"  );
   ImportersToolboxDebug = Journal_Register( Debug_Type, (Name)"ImportersToolboxDebug" );
   ImportersToolboxError = Journal_Register( Error_Type, (Name)"ImportersToolboxError" );

   /* Register your components here.. Components which are not instantied do not need to be registered */
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_GocadProperties_Type,     (Name)"0", _VoxelDataHandler_GocadProperties_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_GocadMaterials_Type,      (Name)"0", _VoxelDataHandler_GocadMaterials_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_Geomodeller_Type,         (Name)"0", _VoxelDataHandler_Geomodeller_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_ASCII_Type,               (Name)"0", _VoxelDataHandler_ASCII_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_GMT_Type,                 (Name)"0", _VoxelDataHandler_GMT_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelDataHandler_VTKStructuredPoints_Type, (Name)"0", _VoxelDataHandler_VTKStructuredPoints_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelParticleLayout_Type,                  (Name)"0", _VoxelParticleLayout_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelFieldVariable_Type,                   (Name)"0", _VoxelFieldVariable_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelFieldVariable_GMT_Type,               (Name)"0", _VoxelFieldVariable_GMT_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, VoxelHeightFieldVariable_Type,             (Name)"0", _VoxelHeightFieldVariable_DefaultNew  );
   Stg_ComponentRegister_Add( componentRegister, SpatialDataFieldVariable_Type,             (Name)"0", _SpatialDataFieldVariable_DefaultNew  );

   /* Register Parents for type checking */
   RegisterParent( VoxelDataHandler_Abstract_Type,            Stg_Component_Type );
   RegisterParent( VoxelDataHandler_GocadAbstract_Type,       VoxelDataHandler_Abstract_Type );
   RegisterParent( VoxelDataHandler_GocadProperties_Type,     VoxelDataHandler_GocadAbstract_Type );
   RegisterParent( VoxelDataHandler_GocadMaterials_Type,      VoxelDataHandler_GocadAbstract_Type );
   RegisterParent( VoxelDataHandler_ASCII_Type,               VoxelDataHandler_Abstract_Type );
   RegisterParent( VoxelDataHandler_GMT_Type,                 VoxelDataHandler_Abstract_Type );
   RegisterParent( VoxelDataHandler_VTKStructuredPoints_Type, VoxelDataHandler_Abstract_Type );
   RegisterParent( VoxelFieldVariable_Type,                   FieldVariable_Type );
   RegisterParent( VoxelFieldVariable_GMT_Type,               VoxelFieldVariable_Type );
   RegisterParent( VoxelParticleLayout_Type,                  GlobalParticleLayout_Type );
   RegisterParent( VoxelDataHandler_Geomodeller_Type,         VoxelDataHandler_Abstract_Type );
   RegisterParent( VoxelHeightFieldVariable_Type,             VoxelFieldVariable_Type );
   RegisterParent( SpatialDataFieldVariable_Type,             FieldVariable_Type );


   return True;
}



