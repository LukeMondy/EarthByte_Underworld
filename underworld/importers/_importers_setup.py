import underworld

def voxelField(filename="", propertyID=""):
	"""
	Creates a voxel field variable from a GoCAD voxel property field.

	Args:
		filename (String): point to the file location
		propertyID (String): name of the voxel property field
	Returns:
		None

	"""
	if filename == "":
		underworld.utils.sendError("Must specify a CSV filename")
	if propertyID == "":
		underworld.utils.sendError("Must specify the property field name")

	globalDict = underworld.dictionary.GetDictionary()

	voxelDatahandler = underworld.dictionary.UpdateDictWithComponent(	globalDict,
																		name = "Materials_Voxel_Datahandler",
																		Type = "VoxelDataHandler_GocadProperties",
																		filename = filename,
																		PropertyName = propertyID,
																		mapIAxisToStgAxis = "X",
																		mapJAxisToStgAxis = "Y",
																		mapKAxisToStgAxis = "Z"
																		)
	voxelField = underworld.dictionary.UpdateDictWithComponent(	globalDict,
																name = "Materials_VoxelField",
																Type = "VoxelFieldVariable",
																VoxelDataHandler = "Materials_Voxel_Datahandler"
																)
	voxelPpc = underworld.dictionary.UpdateDictWithComponent(	globalDict,
																name = "Materials_VoxelFieldPpc",
																Type = "Ppc_Variable",
																FieldVariable = "Materials_VoxelField"
																)
	return [voxelDatahandler, voxelField, voxelPpc]
