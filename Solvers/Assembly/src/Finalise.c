#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "Finalise.h"

#include <stdio.h>

Bool Solvers_Assembly_Finalise( void ) {
	Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	
	Stream_IndentBranch( StgFEM_Debug );
	Stream_UnIndentBranch( StgFEM_Debug );
	return True;
}


