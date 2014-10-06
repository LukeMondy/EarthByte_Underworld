#ifndef __Spherical_Toolbox_h__
#define __Spherical_Toolbox_h__

extern const Type Spherical_Toolbox_Type;

typedef struct {
	__Codelet
} Spherical_Toolbox;

void Spherical_Toolbox_Initialise(PluginsManager* pluginsManager, int* argc, char*** argv);

Index Spherical_Register( PluginsManager* pluginsManager );

void Spherical_Toolbox_Finalise( PluginsManager* pluginsManager );

#endif	
