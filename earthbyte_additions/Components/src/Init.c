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
#include <StgFEM/StgFEM.h>

#include <PICellerator/PICellerator.h>

#include "Components.h"

Stream* earthbyte_additionsInfo  = NULL;
Stream* earthbyte_additionsDebug = NULL;
Stream* earthbyte_additionsError = NULL;

Bool earthbyte_additions_Components_Init()
{
    Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

    Journal_Printf( Journal_Register( Debug_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

    /* === Register components === */
        /* PPC additions */
        Stg_ComponentRegister_Add( componentRegister, Ppc_Material_Condition_Type,     (Name)"0", _Ppc_Material_Condition_DefaultNew  );
        Stg_ComponentRegister_Add( componentRegister, Ppc_Current_Time_Type,           (Name)"0", _Ppc_Current_Time_DefaultNew  );  
        Stg_ComponentRegister_Add( componentRegister, Ppc_Random_Type,                 (Name)"0", _Ppc_Random_DefaultNew  );  
        Stg_ComponentRegister_Add( componentRegister, Ppc_Thermal_Profile_Type,        (Name)"0", _Ppc_Thermal_Profile_DefaultNew  );  
        Stg_ComponentRegister_Add( componentRegister, Ppc_Auto_Thermal_Profile_Type,   (Name)"0", _Ppc_Auto_Thermal_Profile_DefaultNew  );  
        Stg_ComponentRegister_Add( componentRegister, Ppc_Melt_Polynomial_Type,		   (Name)"0", _Ppc_Melt_Polynomial_DefaultNew  );  
		Stg_ComponentRegister_Add( componentRegister, Ppc_PartialMelt_Limited_Type,    (Name)"0", _Ppc_PartialMelt_Limited_DefaultNew  );
        Stg_ComponentRegister_Add( componentRegister, Ppc_Curie_Condition_Type,        (Name)"0", _Ppc_Curie_Condition_DefaultNew  );

        /* Underworld Utils additions */
        Stg_ComponentRegister_Add( componentRegister, PressureTemperatureStrainRateOutput_Type, (Name)"0", _PressureTemperatureStrainRateOutput_DefaultNew  );
   
   
    /* === Register Parents for type checking === */
        /* PPC additions */
        RegisterParent( Ppc_Material_Condition_Type,       Ppc_Type );
        RegisterParent( Ppc_Current_Time_Type,             Ppc_Type );
        RegisterParent( Ppc_Random_Type,                   Ppc_Type );
        RegisterParent( Ppc_Thermal_Profile_Type,          Ppc_Type );
        RegisterParent( Ppc_Auto_Thermal_Profile_Type,     Ppc_Type );
		RegisterParent( Ppc_Melt_Polynomial_Type,  		   Ppc_Type );
        RegisterParent( Ppc_PartialMelt_Limited_Type,      Ppc_Type );
        RegisterParent( Ppc_Curie_Condition_Type,          Ppc_Type );

        /* Underworld Utils additions */
        RegisterParent( PressureTemperatureStrainRateOutput_Type, SwarmOutput_Type );

    return True;
}



