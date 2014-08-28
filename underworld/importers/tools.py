import underworld
import numpy

def LayersFromCSV(CSVfile="", voxel_field="Materials_VoxelField"):
	"""
	1. Imports all the material properties from a CSV file (ignores first row).
		CSV:

		| lithology     | material-index | conductivity (W/mK) | heat production (uW/m^3) |
		-----------------------------------------------------------------------------------
		| 'sediments'   |        1       |         1.75        |        0.70              |
		| 'basement'    |        2       |         3.20        |        2.25              |
		etc...

	2. Constructs a python dictionary with material properties for each lithology.

	3. Sets up shapes and materials corresponding to a voxel field

	Args:
		CSVfile (String): location of CSV file
		voxel_field: name of voxel field in the dictionary (default: "Materials_VoxelField")
	Returns:
		properties (Dict): 'lithology' : [lithology_index, conductivity, heat_production]
		lithology (String): list of names to query material properties from the dictionary ^
	"""
	if CSVfile == "":
		print "You must specify an input CSV file"

	lithology, idx = numpy.loadtxt(CSVfile, delimiter=',', dtype=str, usecols=(0,1), skiprows=1, unpack=True)
	idx, k, HP = numpy.loadtxt(CSVfile, delimiter=',', dtype=float, usecols=(1,2,3), skiprows=1, unpack=True)
	HP = HP*1e-6
	index = 0
	properties = {}
	for ref in lithology:
		properties[ref] = [float(idx[index]), float(k[index]), float(HP[index])]

		shape_dict = underworld.shape.setup.VoxelFieldShape(	shape_name = str(ref)+"_shape",
																LowerLimit = float(idx[index]) -0.5,
																UpperLimit = float(idx[index]) +0.5
																)
		# shape_dict = uw.dictionary.UpdateDictWithComponent(name = str(ref)+"_shape",
		# 	Type = "FieldValueShape",
		# 	ValueField = voxel_field,
		# 	LowerLimit = float(idx[index]) -0.5,
		# 	UpperLimit = float(idx[index]) +0.5
		# 	)
		material_dict = underworld.dictionary.UpdateDictWithComponent(	name = str(ref)+"_material",
																		Type = "Material",
																		Shape = str(ref)+"_shape",
																		thermalConductivity = float(k[index]),
																		heatProduction = float(HP[index])
																		)

		index += 1

	return properties, lithology