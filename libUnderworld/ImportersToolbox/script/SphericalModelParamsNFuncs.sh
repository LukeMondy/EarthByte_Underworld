#!/bin/bash
#
# Underworld Importer Geothermal UW model parameters
#
# Copyright (C) 2012, Monash University
#
# Contact: John Mansour


# Define input variables variables. Add new items as required.
# ------------------------------------------------------------

flag="name"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="A name for your simulation"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTNAME
unset flag; unset item; unset testfntext;

flag="GMTFilename"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Filename for the GMT dataset"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$NODEFAULT 
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if [ \$# -lt 1 ]; then
      echo -e "\nError: You must provide a filename"
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

flag="elementResLat"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements in Latitude"
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTRES
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 1 ]; then
      echo -e "\nError: A simulation element resolution must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="elementResLong"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements in Longitude"
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTRES
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 1 ]; then
      echo -e "\nError: A simulation element resolution must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="elementResR"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements Radially"
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTRES
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9]+$ ]] ; then
      echo -e "\nError: \$1 is invalid.  You must provide a positive integer number."
      return 1
   fi
   if [ \$# -lt 1 ]; then
      echo -e "\nError: A simulation element resolution must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;

flag="minLat"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain Latitude."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="-30"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

flag="maxLat"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain Latitude."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="30"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

flag="minLong"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain Longitude."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="-30"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

flag="maxLong"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain Longitude."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="30"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

flag="minR"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain Radial."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=".5"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

flag="maxR"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain Radial."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="1.0"
testfntext=$(cat <<EOF                                          # define a test function
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
unset flag; unset item; unset testfntext; unset defaultfntext;

# Define required functions.
# --------------------------

function writeuwheaderopen(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local text=$(
cat << 'EOF'
<?xml version="1.0"?>
<!DOCTYPE StGermainData SYSTEM "stgermain.dtd">

<StGermainData xmlns="http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003">
EOF
)
echo ${text} >> $1
unset IFS
}

function writeuwheaderclose(){
echo "</StGermainData>" >> $1
}

function writeinclude(){
echo "  " "<include>" $2 "</include>" >> $1
}

function writeimports(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local text=$(
cat << 'EOF'
   <list name="import" mergeType="merge">
      <param> ImportersToolbox </param>
   </list>
EOF
)
echo ${text} >> $1
unset IFS
}

function writeparameter(){
eval $(printf "%$4s<param name=\"%s\"> %s </param>\n" "" $2 $3 >>$1 )
}
function writecommentedparameter(){
eval $(printf "%$4s<!--param name=\"%s\"> %s </param-->\n" "" $2 $3 >>$1 )
}

function writeparameterformat(){
eval $(printf "%$6s<param name=\"%$4s\"> %$5s </param>\n" "" $2 $3 >>$1 )
}

function writecomment(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local space="                                        "
echo ${space:0:$3-1} "<!--" $2 "-->" >>$1 
unset IFS
}

function writeopenstructreplace(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local space="                                        "
echo ${space:0:$3-1} "<struct name=\"${2}\" mergeType=\"replace\">" >>$1 
unset IFS
}

function writeopenstructmerge(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local space="                                        "
echo ${space:0:$3-1} "<struct name=\"${2}\" mergeType=\"merge\">" >>$1 
unset IFS
}

function writeopenstruct(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local space="                                        "
echo ${space:0:$3-1} "<struct name=\"${2}\">" >>$1 
unset IFS
}

function writeclosestruct(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local space="                                        "
echo ${space:0:$2-1} "</struct>" >>$1 
unset IFS
}



function createGMTVoxelPropertyField() {  # outputFilename providedObjectUniqueName gmtfilename dataPos  indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct $1 "$2_Voxel_Datahandler" $farg
writeparameter $1 "Type" "VoxelDataHandler_GMT" $fargind
writeparameter $1 "filename" $3 $fargind
writeparameter $1 "DataPos" $4 $fargind
writeparameter $1 "DataType" "float" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelField" $farg
writeparameter $1 "Type" "VoxelFieldVariable_GMT" $fargind
writeparameter $1 "VoxelDataHandler" "$2_Voxel_Datahandler" $fargind
writecommentedparameter $1  "PlateDepth" "0.2" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelFieldPpc" $farg
writeparameter $1 "Type" "Ppc_Variable" $fargind
writeparameter $1 "FieldVariable" "$2_VoxelField" $fargind
writeclosestruct $1 $farg 
unset IFS
}




