#!/bin/bash
#
# Underworld GPlates model generator
#
# Copyright (C) 2013, Monash University
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
DEFAULTNAME=GPlates
DEFAULTRES=16
DEFAULTLONG=5
DEFAULTLAT=4
DEFAULTPLATEID=3
DEFAULTEXPAND=0

USAGETOPINFO="Usage: `basename $0` { --name="Foo" --elementResLat=32 ... }"

# Include any required script files
# ---------------------------------

source ${DIR}/BasicParamsNFuncs.sh
source ${DIR}/SphericalModelParamsNFuncs.sh
source ${DIR}/RuntimeParamsNFuncs.sh

# Define input variables variables. Add new items as required.
# ------------------------------------------------------------

flag="latDataColumn"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Column containing your latitude velocity data.  Choose '0' if you do not wish to use velocity data."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTLAT
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 0 ]; then
      echo -e "\nError: A non-negative value must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="longDataColumn"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Column containing your longitude velocity data.  Choose '0' if you do not wish to use velocity data."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTLONG
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 0 ]; then
      echo -e "\nError: A non-negative value must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="plateidDataColumn"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Column containing your plateID data.  Choose '0' if you do not wish to use plateID data."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTPLATEID
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 0 ]; then
      echo -e "\nError: A non-negative value must be provided"
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
writeinclude ${filename} "ImportersToolbox/gplateVoxelGMTGlobeBaseApp.xml"
echo "" >> ${filename}
writecomment ${filename} "User Defined Parameters" 3
writecomment ${filename} "#######################" 3
echo "" >> ${filename}
writecomment ${filename} " radius (> 0) " 3
writeparameter ${filename} "elementResR" ${elementResRArg} 3
writeparameter ${filename} "minR" ${minRArg} 3
writeparameter ${filename} "maxR" ${maxRArg} 3
echo "" >> ${filename}
writecomment ${filename} " longitude range is [-180,180) " 3
writeparameter ${filename} "elementResLong" ${elementResLongArg} 3
writeparameter ${filename} "minLong" ${minLongArg} 3
writeparameter ${filename} "maxLong" ${maxLongArg} 3
echo "" >> ${filename}
writecomment ${filename} " latitude - range is (-90,90) " 3
writeparameter ${filename} "elementResLat" ${elementResLatArg} 3
writeparameter ${filename} "minLat" ${minLatArg} 3
writeparameter ${filename} "maxLat" ${maxLatArg} 3
echo "" >> ${filename}

writecomment ${filename} " Gravity " 3
writeparameter ${filename} "gravity" "10" 3
echo "" >> ${filename}

writecomment ${filename} " Simulation control " 3
writeparameter ${filename} "maxTimeSteps" "0" 3
writeparameter ${filename} "dumpEvery" "1" 3
writeparameter ${filename} "checkpointEvery" "1" 3
writeparameter ${filename} "outputPath" "./output" 3
echo "" >> ${filename}


writecomment ${filename} "Boundary Conditions" 3
if [ $longDataColumnArg -ne 0 ]; then
   writeinclude ${filename} "ImportersToolbox/gplateGlobe_BCs.xml"
else
   writeinclude ${filename} "gplateGlobe_BCFreeSlip.xml"
fi
echo "" >> ${filename}

writeopenstructmerge ${filename} "components" 3

if [ $longDataColumnArg -ne 0 ]; then
   createGMTVoxelPropertyField ${filename} "Longitude" ${GMTFilenameArg} $longDataColumnArg 6
   echo "" >> ${filename}
fi
if [ $latDataColumnArg -ne 0 ]; then
   createGMTVoxelPropertyField ${filename} "Latitude" ${GMTFilenameArg} $latDataColumnArg 6
   echo "" >> ${filename}
fi
if [ $plateidDataColumnArg -ne 0 ]; then
   createGMTVoxelPropertyField ${filename} "PlateID" ${GMTFilenameArg} $plateidDataColumnArg 6
   echo "" >> ${filename}
   writecomment ${filename} "Add Plate Materials" 6
   $DIR/getShapesFromVelocityFile.pl < ${GMTFilenameArg} >>  ${filename}
fi

writeclosestruct ${filename} 3

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
