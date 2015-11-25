
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include "Spherical.h"

#include <stdio.h>

/** Initialises this package, then any init for this package
such as streams etc */

Bool Spherical_Init( int* argc, char** argv[] ) {
	/* This init function tells StGermain of all the component types, etc this module contributes. Because it can be linked at compile
	   time or linked in by a toolbox at runtime, we need to make sure it isn't run twice (compiled in and loaded through a toolbox.*/
	if( !ToolboxesManager_IsInitialised( stgToolboxesManager, "Spherical" ) ) {
		int tmp;
		char* directory;

		/*Spherical_Base_Init(argc, argv); */

		Journal_Printf( Journal_Register( DebugStream_Type, (Name)"Context"  ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
		tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context" )  );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), 0 );
		Journal_Printf( /* DO NOT CHANGE OR REMOVE */
			Journal_Register( InfoStream_Type, (Name)"Context"  ), 
			"Spherical (Bleeding Edge Geodynamics framework). Copyright (C) 2005 Monash University.\n" );
		Stream_Flush( Journal_Register( InfoStream_Type, (Name)"Context" )  );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, (Name)"Context"  ), tmp );

		/* Add the Spherical path to the global xml path dictionary */
		directory = Memory_Alloc_Array( char, 200, "xmlDirectory" ) ;
		sprintf(directory, "%s%s", LIB_DIR, "/StGermain" );
		XML_IO_Handler_AddDirectory( "Spherical", directory );
		Memory_Free(directory);

		/* Add the plugin path to the global plugin list */
		ModulesManager_AddDirectory( "Spherical", LIB_DIR );
	
		return True;
	}
	return False;
}
