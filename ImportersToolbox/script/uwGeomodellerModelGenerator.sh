#!/bin/bash
#
# Underworld Geomodeller Geothermal model builder
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
DEFAULTNAME=GeomodellerGeothermal
DEFAULTRES=16
DEFAULTTEMPLOWER=300
DEFAULTTEMPUPPER=0
DEFAULTFLUXLOWER=0.0125
DEFAULTCONDUCTIVITY=1.
DEFAULTHEATPRODUCTION=0.
DEFAULTEXPAND=0

USAGETOPINFO="Usage: `basename $0` { --name="Foo" --elementResI=32 ... }"
# Include any required script files
# ---------------------------------

source ${DIR}/BasicParamsNFuncs.sh
source ${DIR}/ModelParamsNFuncs.sh
source ${DIR}/RuntimeParamsNFuncs.sh

# Define input variables variables. Add new items as required.
# ------------------------------------------------------------


flag="defaultConductivity"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Default thermal conductivity to be used in regions where voxets do not specify a value."
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTCONDUCTIVITY
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9.-]+([.][0-9]+)?$ ]] ; then
      echo -e "\nError: \$1 is not a number"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="defaultHeatProduction"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Default heat production to be used in regions where voxets do not specify a value."
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTHEATPRODUCTION
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9.-]+([.][0-9]+)?$ ]] ; then
      echo -e "\nError: \$1 is not a number"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="keyFilename"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of your materials key CSV file.  Data must be in columns:  Name (as per voxet), Conductivity, HeatProduction."
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a materials key filename"
      return 1
   fi
   if [ ! -f \$1 ]; then
      echo -e "\nError: File '\$1' not found"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

function KeyValidIndexTest () {
if ! [[ "$1" =~ ^[0-9,e,E,.]+$ ]] ; then
   return 1
fi
return 0
}

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
   echo "Running Setup Wizard.  If available, default value used for values not provided."
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

# generate arrays with material information from key file
#while IFS=, read var1 var2 var3 || [[ -n "$LINE" ]]; do
# need to use below process to ensure last line will be read!  what a pain.
while IFS=$'\n' read -r LINE || [[ -n "$LINE" ]]; do
   var1=$(echo $LINE | awk -F, '{print $1}')
   var2=$(echo $LINE | awk -F, '{print $2}')
   var3=$(echo $LINE | awk -F, '{print $3}')
   KeyValidIndexTest $var2
   returnVal=$?
   if [ $returnVal -eq 0 ]; then
      materialName=`echo $var1 | sed 's/[- ]/_/g'` # replace spaces and dashes in names with underscores
      lenstring=${#materialName}                   # get max name length
      matNames[${#matNames[*]}]=$materialName 
      matConduct[${#matConduct[*]}]=$var2
      matHeat[${#matHeat[*]}]=$var3
   fi
done <  ${keyFilenameArg}


numMats=${#matNames[*]}
maxLenName=`getMaxVarLength "${matNames[@]}"`
maxLenProps=`getMaxVarLength "${matConduct[@]}" "${matHeat[@]}"`
# Generate model xml

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

writecomment ${filename} "Write material properties" 3
writeparameter ${filename} "defaultConductivity" ${defaultConductivityArg} 3
writeparameter ${filename} "defaultHeatProduction" ${defaultHeatProductionArg} 3
echo "" >> ${filename}
for (( matInd = 0 ; matInd < $numMats ; matInd++ ))
do
   writeparameterformat ${filename} ${matNames[matInd]}_k ${matConduct[matInd]} $[$maxLenName+2] $maxLenProps 3
   writeparameterformat ${filename} ${matNames[matInd]}_q ${matHeat[matInd]} $[$maxLenName+2] $maxLenProps 3
done
echo "" >> ${filename}

writecomment ${filename} "Voxet Files/Props to Load" 3
matvoxfile="MatVoxelFile"
writeparameter ${filename} ${matvoxfile} ${voxetFilenameArg} 3
echo "" >> ${filename}

writecomment ${filename} "Model Construction (do not modify)" 3
writecomment ${filename} "##################################" 3
writeopenstructmerge ${filename} "components" 3
writeopenstruct ${filename} "EverywhereShape" 6
writeparameter ${filename} "Type" "Everywhere" 9
writeclosestruct ${filename} 6
writeopenstruct ${filename} "EverywhereMaterial" 6
writeparameter ${filename} "Type" "Material" 9
writeparameter ${filename} "Shape" "EverywhereShape" 9
writeparameter ${filename} "thermalConductivity" "@defaultConductivity" 9
writeparameter ${filename} "heatProduction" "@defaultHeatProduction" 9
writeclosestruct ${filename} 6
echo "" >> ${filename}
createGeomodellerVoxelPropertyField $filename "Materials" "@${matvoxfile}" "@${matvoxprop}"  6
echo "" >> ${filename}
writecomment ${filename} "Write materials" 6
for (( matInd = 0 ; matInd < $numMats ; matInd++ ))
do
   createShapedMaterial ${filename} "${matNames[matInd]}" "Materials_VoxelField" "@${matNames[matInd]}_k" "@${matNames[matInd]}_q" 6
done
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
