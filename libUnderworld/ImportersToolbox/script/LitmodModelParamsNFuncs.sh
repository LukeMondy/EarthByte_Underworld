#!/bin/bash
#
# Underworld Importer Litmod UW model parameters
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

flag="elementResI"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements along X axis"
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

flag="elementResJ"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements along Y axis"
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

flag="elementResK"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of FEM elements along Z axis"
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

flag="expandDomain"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Set a domain expansion/contraction factor"
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTEXPAND
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


flag="minX"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain X coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$1}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

flag="maxX"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain X coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$2}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

flag="minY"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain Y coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$3}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

flag="maxY"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain Y coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$4}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

flag="minZ"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Minimum domain Z coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$5}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

flag="maxZ"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Maximum domain Z coordinate. If available, default is obtained from vtk file."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVTK"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$VTKFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$6}' 2>/dev/null)
      defaultvalue[$item]=\$__result
   }
EOF
)
eval "$defaultfntext"
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

