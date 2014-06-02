#!/bin/bash
#
# Underworld Litmod model generator
#
# Copyright (C) 2014, Monash University
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
DEFAULTNAME=Litmod
DEFAULTRES=16
DEFAULTEXPAND=1

USAGETOPINFO="Usage: `basename $0` { --name="Foo" --elementResLat=32 ... }"

# Include any required script files
# ---------------------------------

source ${DIR}/BasicParamsNFuncs.sh
source ${DIR}/LitmodModelParamsNFuncs.sh
source ${DIR}/RuntimeParamsNFuncs.sh

# add some parameters

flag="VTKFilename"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Litmod VTK filename"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a vtk filename"
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

flag="DensityDatasetName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="VTK dataset name corresponding to denity field"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a dataset name"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="TemperatureDatasetName"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="VTK dataset name corresponding to temperature field"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a dataset name"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
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
writeinclude ${filename} "ImportersToolbox/LitmodBase3D.xml"
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

writeparameter ${filename} "TemperatureVTKFilename" ${VTKFilenameArg} 3
writeparameter ${filename} "TemperatureVTKDataset"  ${TemperatureDatasetNameArg} 3
writeparameter ${filename} "DensityVTKFilename"     ${VTKFilenameArg} 3
writeparameter ${filename} "DensityVTKDataset"      ${DensityDatasetNameArg} 3
echo "" >> ${filename}

writeparameter ${filename} "gravity" "10" 3
writeparameter ${filename} "eta0" "1" 3
writeparameter ${filename} "cohesion" "1.0e10" 3
writeparameter ${filename} "frictionCoefficient" "0.1" 3
echo "" >> ${filename}

writecomment ${filename} " Simulation control " 3
writeparameter ${filename} "maxTimeSteps" "0" 3
writeparameter ${filename} "dumpEvery" "1" 3
writeparameter ${filename} "checkpointEvery" "1" 3
writeparameter ${filename} "outputPath" "./output" 3
echo "" >> ${filename}


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
