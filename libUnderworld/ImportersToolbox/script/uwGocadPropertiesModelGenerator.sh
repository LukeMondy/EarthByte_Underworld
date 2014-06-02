#!/bin/bash
#
# Underworld Gocad Geothermal model builder
#
# Copyright (C) 2012, Monash University
#
# Contact: John Mansour

# Get actual script location, including dereferenced symbolic links.
# Note that this script relies on other scripts found in its original location.
# As such, you may not move this script, but instead provide a symbolic link if required. 

SOURCE="${BASH_SOURCE[0]}"
DIR="$( dirname "$SOURCE" )"
while [ -h "$SOURCE" ]
do 
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
  DIR="$( cd -P "$( dirname "$SOURCE"  )" && pwd )"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# set a few defaults
NODEFAULT=
DEFAULTNAME=GocadGeothermal
DEFAULTRES=16
DEFAULTTEMPLOWER=300
DEFAULTTEMPUPPER=0
DEFAULTFLUXLOWER=0.0125
DEFAULTEXPAND=0
USAGETOPINFO="Usage: `basename $0` { --name="Foo" --elementResI=32 ... }"

# Include any required script files
# ---------------------------------

source ${DIR}/BasicParamsNFuncs.sh
source ${DIR}/ModelParamsNFuncs.sh
source ${DIR}/RuntimeParamsNFuncs.sh

# Define input variables variables. Add new items as required.
# ------------------------------------------------------------

flag="conductivityPropName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of voxet property corresponding to the conductivity"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a conductivity voxet property name"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;


flag="heatPropName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of voxet property corresponding to the heat generation"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
unset flag; unset item; unset testfntext;


flag="densityPropName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of voxet property corresponding to density (not generally required)"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]="NOTSET"
unset flag; unset item; unset testfntext;

# Run
# =======================================================================
#

# Ensure there's enough command line arguments...
EXPECTED_ARGS=1
E_BADARGS=65
if [ $# -lt $EXPECTED_ARGS ]; then
   usage
   exit $E_BADARGS
fi

# first check for all argumentless options
evalNoOptArgs $@

# run help command now if requested
if [ $helpSet ]; then
   usage
   exit
fi

# if running wizard, print info
if [ $wizardSet ]; then
   echo ""
   echo "Running Setup Wizard.  If available, default value used for values not entered."
   echo ""
fi

# now check for all argumented options.  use wizard if requested.
allargedoptions=("${takesargument[@]}" "${uwtakesargument[@]}")
for option in "${allargedoptions[@]}";
do
   readarg $option $@
done

usedValues


filename=${nameArg}.xml
echo ""
echo "Writing Underworld model file" ${filename}

if [ -f "${filename}" ]
then
   mv ${filename} ${filename}.orig
fi

writeuwheaderopen ${filename} 
echo "" >> ${filename}
writeimports ${filename}
echo "" >> ${filename}
writeinclude ${filename} "Geothermal/HeatFlowApp.xml"
writeinclude ${filename} "Geothermal/_HeatFlowSolver_Ppc_UW.xml"
echo "" >> ${filename}
writecomment ${filename} "User Defined Parameters" 3
writecomment ${filename} "#######################" 3
echo "" >> ${filename}

writecomment ${filename} "Mesh Configuration" 3
writeparameter ${filename} "dim" 3 3
echo "" >> ${filename}
writeparameter ${filename} "elementResI" ${elementResIArg} 3
writeparameter ${filename} "elementResJ" ${elementResJArg} 3
writeparameter ${filename} "elementResK" ${elementResKArg} 3
echo "" >> ${filename}
writecomment ${filename} "Domain" 3
writeparameter ${filename} "minX" ${minXArg} 3
writeparameter ${filename} "minY" ${minYArg} 3
writeparameter ${filename} "minZ" ${minZArg} 3
writeparameter ${filename} "maxX" ${maxXArg} 3
writeparameter ${filename} "maxY" ${maxYArg} 3
writeparameter ${filename} "maxZ" ${maxZArg} 3
echo "" >> ${filename}

writecomment ${filename} "Boundary Conditions" 3
writeparameter ${filename} "upperTemp" ${upperBoundaryTempArg} 3
if [ $useLowerFluxBCSet ]; then
   writeparameter ${filename} "lowerFlux" ${lowerBoundaryFluxArg} 3
else
   writeparameter ${filename} "lowerTemp" ${lowerBoundaryTempArg} 3
fi
echo "" >> ${filename}

writecomment ${filename} "Data Output Path" 3
writeparameter ${filename} "outputPath" "./output_"${nameArg} 3
echo "" >> ${filename}
writecomment ${filename} "Model Construction (do not modify)" 3
writecomment ${filename} "##################################" 3
writeopenstructmerge ${filename} "components" 3
echo "" >> ${filename}

createGocadVoxelPropertyField ${filename} "Conductivity" ${voxetFilenameArg} ${conductivityPropNameArg} 6

heatFieldToUse="0."
if [ $heatPropNameArg ]; then
   heatFieldToUse="Heat_VoxelFieldPpc"
   createGocadVoxelPropertyField ${filename} "Heat" ${voxetFilenameArg} ${heatPropNameArg} 6
fi

if [ ${densityPropNameArg} != "NOTSET"  ]; then
   createGocadVoxelPropertyField ${filename} "Density" ${voxetFilenameArg} ${densityPropNameArg} 6
fi

writeopenstruct ${filename} "EverywhereShape" 6
writeparameter ${filename} "Type" "Everywhere" 9
writeclosestruct ${filename} 6
writeopenstruct ${filename} "EverywhereMaterial" 6
writeparameter ${filename} "Type" "Material" 9
writeparameter ${filename} "Shape" "EverywhereShape" 9
writeparameter ${filename} "thermalConductivity" "Conductivity_VoxelFieldPpc" 9
writeparameter ${filename} "heatProduction" ${heatFieldToUse} 9
writeclosestruct ${filename} 6
echo "" >> ${filename}
if [ $useLowerFluxBCSet ]; then
  createFluxTerm ${filename} 6
fi 
writeclosestruct ${filename} 3
echo "" >> ${filename}


if [ $useLowerFluxBCSet ]; then
   writeUpperOnlyTemperatureBC ${filename}
else
   writeUpperAndLowerTemperatureBC ${filename}
fi 

writeuwheaderclose ${filename}


echo "Done writing model file"

# Run model 

if [ $runSet ]; then
   echo "Now running model with nproc =" $nprocsArg
   if [ $nprocsArg -gt 1 ]; then
      CMD="${mpirunArg} -n ${nprocsArg} ${binArg} ${filename}"
   else
      CMD="${binArg} ${filename}"
   fi
   echo "Command used:"
   echo $CMD
   $CMD
fi
