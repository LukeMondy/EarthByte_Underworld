#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Spherical.h"

#include <stdio.h>

Bool Spherical_Finalise( void ) {
	if( ToolboxesManager_IsInitialised( stgToolboxesManager, "Spherical" ) ) {
		Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

		/*Spherical_Base_Finalise();*/

		return True;
	} else {
		return False;
	}
}


