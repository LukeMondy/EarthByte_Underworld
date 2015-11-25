#!/bin/bash
#
# Underworld Gocad Temperature dependent Geothermal model builder
# 
# Conductivity defined as:
#
#   k = ((T0*TCrit)/(TCrit-T0))*(K0-KCrit)*(1/temperature - 1/TCrit) + KCrit;
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
DEFAULTNAME=GocadTempDepCondGeothermal
DEFAULTRES=16
DEFAULTTEMPLOWER=600
DEFAULTTEMPUPPER=300
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

flag="materialsPropName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of voxet property corresponding to the material definitions"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a materials voxet property name"
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
description[$item]="Name of your materials key CSV file.  Data be in format:  Name, Index, T0, TCrit, K0, KCrit, HeatProduction.  First line a header."
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
if ! [[ "$1" =~ ^[0-9]+$ ]] ; then
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
#while IFS=, read var1 var2 var3 var4 var5 var6 var7
# need to use below process to ensure last line will be read!  what a pain.
while IFS=$'\n' read -r LINE || [[ -n "$LINE" ]]; do
   var1=$(echo $LINE | awk -F, '{print $1}')
   var2=$(echo $LINE | awk -F, '{print $2}')
   var3=$(echo $LINE | awk -F, '{print $3}')
   var4=$(echo $LINE | awk -F, '{print $4}')
   var5=$(echo $LINE | awk -F, '{print $5}')
   var6=$(echo $LINE | awk -F, '{print $6}')
   var7=$(echo $LINE | awk -F, '{print $7}')
   KeyValidIndexTest $var2
   returnVal=$?
   if [ $returnVal -eq 0 ]; then
      materialName=`echo $var1 | sed 's/[- ]/_/g'` # replace spaces and dashes in names with underscores
      lenstring=${#materialName}                   # get max name length
      matNames[${#matNames[*]}]=$materialName 
      matIndex[${#matIndex[*]}]=$var2 
      matConduct_t0[${#matConduct_t0[*]}]=$var3 
      matConduct_tC[${#matConduct_tC[*]}]=$var4 
      matConduct_k0[${#matConduct_k0[*]}]=$var5 
      matConduct_kC[${#matConduct_kC[*]}]=$var6 
      matHeat[${#matHeat[*]}]=`echo $var7 | sed 's/[^0-9\.\-]//g'`
   fi 
done <  ${keyFilenameArg}

numMats=${#matNames[*]}
maxLenName=`getMaxVarLength "${matNames[@]}"`
maxLenProps=`getMaxVarLength "${matConduct_t0[@]}" "${matConduct_tC[@]}" "${matConduct_k0[@]}" "${matConduct_kC[@]}" "${matHeat[@]}"`
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
   writeparameterformat ${filename} ${matNames[matInd]}_t0 ${matConduct_t0[matInd]} $[$maxLenName+2] $maxLenProps 3
   writeparameterformat ${filename} ${matNames[matInd]}_tC ${matConduct_tC[matInd]} $[$maxLenName+2] $maxLenProps 3
   writeparameterformat ${filename} ${matNames[matInd]}_k0 ${matConduct_k0[matInd]} $[$maxLenName+2] $maxLenProps 3
   writeparameterformat ${filename} ${matNames[matInd]}_kC ${matConduct_kC[matInd]} $[$maxLenName+2] $maxLenProps 3
   writeparameterformat ${filename} ${matNames[matInd]}_q ${matHeat[matInd]} $[$maxLenName+2] $maxLenProps 3
   echo "" >> ${filename}
done
echo "" >> ${filename}

writecomment ${filename} "Voxet Files/Props to Load" 3
matvoxfile="MatVoxelFile"
matvoxprop="MatVoxelProperty"
writeparameter ${filename} ${matvoxfile} ${voxetFilenameArg} 3
writeparameter ${filename} ${matvoxprop} ${materialsPropNameArg} 3
echo "" >> ${filename}

writecomment ${filename} "Model Construction (do not modify)" 3
writecomment ${filename} "##################################" 3
writeopenstructmerge ${filename} "components" 3
writeopenstructmerge ${filename} "energySLE" 6
writeparameter ${filename} "isNonLinear" "yes" 9
writeclosestruct ${filename} 6
writeopenstruct ${filename} "EverywhereShape" 6
writeparameter ${filename} "Type" "Everywhere" 9
writeclosestruct ${filename} 6
writeopenstruct ${filename} "EverywhereMaterial" 6
writeparameter ${filename} "Type" "Material" 9
writeparameter ${filename} "Shape" "EverywhereShape" 9
writeparameter ${filename} "thermalConductivity" "2." 9
writeparameter ${filename} "heatProduction" "0." 9
writeclosestruct ${filename} 6
echo "" >> ${filename}
createGocadVoxelPropertyField $filename "Materials" "@${matvoxfile}" "@${matvoxprop}"  6
echo "" >> ${filename}
writecomment ${filename} "Write materials" 6
for (( matInd = 0 ; matInd < $numMats ; matInd++ ))
do
   createTempDepAltPpc      ${filename} "${matNames[matInd]}_k" "@${matNames[matInd]}_t0" "@${matNames[matInd]}_tC" "@${matNames[matInd]}_k0" "@${matNames[matInd]}_kC" 6
   createFieldShapeMaterial ${filename} "${matNames[matInd]}_Shape" "Materials_VoxelField" ${matIndex[matInd]} "${matNames[matInd]}_k" "@${matNames[matInd]}_q" 6
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

writeTempIC ${filename}

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
