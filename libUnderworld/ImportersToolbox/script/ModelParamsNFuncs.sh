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

flag="voxetFilename"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Filename for the voxet dataset"
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

flag="lowerBoundaryTemp"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Lower boundary (minZ) temperature BC value. This value is ignored if --useLowerFluxBC is set."
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTTEMPLOWER
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

flag="upperBoundaryTemp"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Upper boundary (maxZ) temperature BC value"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTTEMPUPPER
testfntext=$(cat <<EOF                                       # define a test function
   function ${flag}Test () {
   if ! [[ "\$1" =~ ^[0-9.-]+([.][0-9]+)?$ ]] ; then
      echo -e "\nError: \$1 is not a number"
      return 1
   fi
   if [ \$# -lt 1 ]; then
      echo -e "\nError: An upper boundary temperature must be provided"
      return 1
   fi
   return 0
   }
EOF
)
eval "$testfntext"
unset flag; unset item; unset testfntext;


flag="useLowerFluxBC"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Use a lower boundary heat flux (Neumann) boundary condition?"
takesnoargument[${#takesnoargument[*]}]=$item
unset flag; unset item; unset testfntext;


flag="lowerBoundaryFlux"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Lower boundary (minZ) heat flux value.  This value is ignored if --useLowerFluxBC is not set."
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=$DEFAULTFLUXLOWER
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
description[$item]="Minimum domain X coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$1}' 2>/dev/null)
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
description[$item]="Maximum domain X coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$2}' 2>/dev/null)
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
description[$item]="Minimum domain Y coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$3}' 2>/dev/null)
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
description[$item]="Maximum domain Y coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$4}' 2>/dev/null)
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
description[$item]="Minimum domain Z coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$5}' 2>/dev/null)
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
description[$item]="Maximum domain Z coordinate. If available, default is value from voxet definition."
uwtakesargument[${#uwtakesargument[*]}]=$item
defaultvalue[$item]="FromVoxet"
defaultfntext=$(cat <<EOF                                       # define a function to get a default
   function ${flag}Default () {
      local __result=\$($DIR/uwGetDomainRange.py --filename=\${argoptionsset[\$voxetFilename]} --domFactor=\${argoptionsset[\$expandDomain]} 2>/dev/null| awk '{print \$6}' 2>/dev/null)
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

                                            #      1                    2                   3           4             5 
function createGocadVoxelPropertyField() {  # outputFilename providedObjectUniqueName voxetfilename propertyName  indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct $1 "$2_Voxel_Datahandler" $farg
writeparameter $1 "Type" "VoxelDataHandler_GocadProperties" $fargind
writeparameter $1 "filename" $3 $fargind
writeparameter $1 "PropertyName" $4 $fargind
writeparameter $1 "mapIAxisToStgAxis" "X" $fargind
writeparameter $1 "mapJAxisToStgAxis" "Y" $fargind
writeparameter $1 "mapKAxisToStgAxis" "Z" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelField" $farg
writeparameter $1 "Type" "VoxelFieldVariable" $fargind
writeparameter $1 "VoxelDataHandler" "$2_Voxel_Datahandler" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelFieldPpc" $farg
writeparameter $1 "Type" "Ppc_Variable" $fargind
writeparameter $1 "FieldVariable" "$2_VoxelField" $fargind
writeclosestruct $1 $farg 
unset IFS
}

                                                #      1                    2                   3           4             5 
function createGeomodellerVoxelPropertyField() {  # outputFilename providedObjectUniqueName voxetfilename propertyName  indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct $1 "$2_Voxel_Datahandler" $farg
writeparameter $1 "Type" "VoxelDataHandler_Geomodeller" $fargind
writeparameter $1 "filename" $3 $fargind
writeparameter $1 "mapIAxisToStgAxis" "X" $fargind
writeparameter $1 "mapJAxisToStgAxis" "Y" $fargind
writeparameter $1 "mapKAxisToStgAxis" "Z" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelField" $farg
writeparameter $1 "Type" "VoxelFieldVariable" $fargind
writeparameter $1 "VoxelDataHandler" "$2_Voxel_Datahandler" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "$2_VoxelFieldPpc" $farg
writeparameter $1 "Type" "Ppc_Variable" $fargind
writeparameter $1 "FieldVariable" "$2_VoxelField" $fargind
writeclosestruct $1 $farg 
unset IFS
}

                                       #      1                     2                3           4            5        6       7
function createFieldShapeMaterial() {  # outputFilename providedObjectUniqueName fieldName materialIndex matConduct matHeat indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct  $1 "$2_FieldValueShape" $farg
writeparameter   $1 "Type" "FieldValueShape" $fargind
writeparameter   $1 "ValueField" $3 $fargind
writeparameter   $1 "LowerLimit" $[$4-1].5 $fargind
writeparameter   $1 "UpperLimit" $[$4].5 $fargind
writeclosestruct $1 $farg
createThermalMaterial $1 $2 $2_FieldValueShape $5 $6 $7
unset IFS
}

                                   #      1                     2                3          4        5       7
function createShapedMaterial() {  # outputFilename providedObjectUniqueName fieldName matConduct matHeat indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct  $1 "$2" $farg
writeparameter   $1 "Type" "ShapedMaterial" $fargind
writeparameter   $1 "MaterialIndexField" $3 $fargind
writeparameter   $1 "thermalConductivity" $4 $fargind
writeparameter   $1 "heatProduction" $5 $fargind
writeclosestruct $1 $farg
unset IFS
}

                                    #      1                     2                  3           4         5         6
function createThermalMaterial() {  # outputFilename providedObjectUniqueName materialShape matConduct matHeat indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct  $1 "$2_Material" $farg
writeparameter   $1 "Type" "Material" $fargind
writeparameter   $1 "Shape" $3 $fargind
writeparameter   $1 "thermalConductivity" $4 $fargind
writeparameter   $1 "heatProduction" $5 $fargind
writeclosestruct $1 $farg
}

function createTempDepAltPpc() {  # outputFilename providedObjectUniqueName t0 tC k0 kC
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct  $1 "$2" $farg
writeparameter   $1 "Type" "Ppc_TempDependentDiffusivityAlt" $fargind
writeparameter   $1 "T0"    $3 $fargind
writeparameter   $1 "TCrit" $4 $fargind
writeparameter   $1 "K0"    $5 $fargind
writeparameter   $1 "KCrit" $6 $fargind
writeclosestruct $1 $farg
}

function writeUpperOnlyTemperatureBC(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local text=$(
cat << 'EOF'
   <struct name="temperatureBCs" mergeType="replace">
      <param name="type">CompositeVC</param>
      <list name="vcList">
         <struct>
            <param name="type">WallVC</param>
            <param name="wall">MaxK</param>
            <list name="variables">
               <struct>
                  <param name="name">temperature</param>
                  <param name="type">double</param>
                  <param name="value">@upperTemp</param>
               </struct>
            </list>
         </struct>
      </list>
   </struct>
EOF
)
echo ${text} >> $1
unset IFS
}


function writeUpperAndLowerTemperatureBC(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local text=$(
cat << 'EOF'
   <struct name="temperatureBCs" mergeType="replace">
      <param name="type">CompositeVC</param>
      <list name="vcList">
         <struct>
            <param name="type">WallVC</param>
            <param name="wall">MinK</param>
            <list name="variables">
               <struct>
                  <param name="name">temperature</param>
                  <param name="type">double</param>
                  <param name="value">@lowerTemp</param>
               </struct>
            </list>
         </struct>
         <struct>
            <param name="type">WallVC</param>
            <param name="wall">MaxK</param>
            <list name="variables">
               <struct>
                  <param name="name">temperature</param>
                  <param name="type">double</param>
                  <param name="value">@upperTemp</param>
               </struct>
            </list>
         </struct>
      </list>
   </struct>
EOF
)
echo ${text} >> $1
unset IFS
}

function writeTempIC(){
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local text=$(
cat << 'EOF'
   <param name="LinearTempIC_TopLayerCoord">@maxZ</param>
   <param name="LinearTempIC_TopLayerBC"> @upperTemp </param>
   <param name="LinearTempIC_BottomLayerCoord">@minZ</param>
   <param name="LinearTempIC_BottomLayerBC"> @lowerTemp </param>
   <param name="VerticalAxis"> 2 </param>
   <struct name="temperatureICs" mergeType="replace">
      <param name="type">CompositeVC</param>
      <list name="vcList">
         <struct>
            <param name="type">AllNodesVC</param>
            <list name="variables">
               <struct>
                  <param name="name">temperature</param>
                  <param name="type">func</param>
                  <param name="value">TemperatureLinear</param>
               </struct>
            </list>
         </struct>
      </list>
   </struct>
EOF
)
echo ${text} >> $1
unset IFS
}
                             #      1             2
function createFluxTerm() {  # outputFilename indentLevel
IFS='%'  # this is required because bash automatically strips whitespace otherwise
local farg=${@: -1}  #final arg, always indent
local fargind=$[$farg+3]
writeopenstruct $1 "FluxTerm" $farg
writeparameter $1 "Type" "VectorSurfaceAssemblyTerm_NA__Fi__ni" $fargind
writeparameter $1 "ForceVector" "fVector" $fargind
writeparameter $1 "functionLabel" "fluxVector" $fargind
writeparameter $1 "Swarm" "gaussBorderSwarm" $fargind
writeparameter $1 "Surface" "MinK" $fargind
writeclosestruct $1 $farg
writeopenstruct $1 "fluxVector" $farg
writeparameter $1 "Type" "Ppc_a_Vector" $fargind
writeparameter $1 "Alpha" "@lowerFlux" $fargind
writeparameter $1 "vi" "0" $fargind
writeparameter $1 "vj" "0" $fargind
writeparameter $1 "vk" "1" $fargind
writeclosestruct $1 $farg
unset IFS
}

