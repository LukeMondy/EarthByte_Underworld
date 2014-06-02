#!/bin/bash
#
# Underworld Importer Geothermal Basic Parameters and Functions
#
# Copyright (C) 2012, Monash University
#
# Contact: John Mansour


# Define minimal parameters for script 
# ------------------------------------------------------------

flag="help"
options[${#options[*]}]=$flag                                 # add item to options list
item=$((${#options[*]}-1))                                    # get index of item
eval ${options[$item]}=$item                                  # declare variable help=0
description[$item]="Print this page"                          # define description
takesnoargument[${#takesnoargument[*]}]=$item                 # add to takesnoargument list
unset flag; unset item; unset testfntext;

flag="wizard"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Run wizard to input your model options"
takesnoargument[${#takesnoargument[*]}]=$item
unset flag; unset item; unset testfntext;


# Define required functions.
# --------------------------

function readarg() {
local pass=0
local providedArgument
currentOption=${options[$1]}
shift
while [ $pass -lt 1 ]; do
   if type ${currentOption}Default &>/dev/null; then            # if default option function provided, first attempt it
      ${currentOption}Default
   fi
   if findoption ${currentOption} $@; then                      # if provided at commandline, use
      providedArgument=`getoptionarg ${currentOption} $@`
   elif [ $wizardSet ]; then                                    # if wizard, use
      echo ${description[$option]}
      if [ "x"${defaultvalue[$option]} = "x" ]; then
         echo -n ": "
      else
         echo -n "("${defaultvalue[$option]}"): "
      fi
      read providedArgument
      echo ""
   fi

   if [ "x"${providedArgument} = "x" ]; then                    # if still nothing set, use default
      if [ ! "x"${defaultvalue[$option]} = "x" ]; then          # if available
         providedArgument=${defaultvalue[$option]}
      fi
   fi
   pass=1                                                       # if no test, default to pass
   if type ${currentOption}Test &>/dev/null; then               # run test if defined
      ${currentOption}Test $providedArgument
      returnVal=$?
      if [ $returnVal -gt 0 ]; then
         pass=0
         if [ ! $wizardSet ]; then
            echo "Set the commandline argument --${currentOption}=, or use the wizard (--wizard)"
            exit 1
         fi
      else
         pass=1
      fi
   fi
   eval ${currentOption}Arg=$providedArgument                   # set argument to for eg elementResIArg=16
   eval ${currentOption}Set=1                                   # set this flag, eg elementResISet=1
   eval argoptionsset[$currentOption]=${providedArgument}       # add to argument list
done
}

function usage() {
IFS='%'  # this is required because bash automatically strips whitespace otherwise
INDENT=40
cat << EOF
$USAGETOPINFO
Options: [ --flag, --flag=, --flag=(default) ]
EOF
for index in ${!takesnoargument[*]}
do
   itemindex=${takesnoargument[$index]}
   FLAG="    --"${options[itemindex]}"                                                                   "
   FULLLINE=${FLAG:0:$INDENT}${description[itemindex]}
   echo $FULLLINE
done
for index in ${!takesargument[*]}
do
   itemindex=${takesargument[$index]}
   if [ "x"${defaultvalue[itemindex]} = "x" ]; then
      FLAG="    --"${options[itemindex]}"=                                                               "
   else
      FLAG="    --"${options[itemindex]}"=("${defaultvalue[itemindex]}")                                 "
   fi
   FULLLINE=${FLAG:0:$INDENT}${description[itemindex]}
   echo $FULLLINE
done
#echo "(underworld and script compatible options)"
for index in ${!uwtakesargument[*]}
do
   itemindex=${uwtakesargument[$index]}
   if [ "x"${defaultvalue[itemindex]} = "x" ]; then
      FLAG="    --"${options[itemindex]}"=                                                               "
   else
      FLAG="    --"${options[itemindex]}"=("${defaultvalue[itemindex]}")                                 "
   fi
   FULLLINE=${FLAG:0:$INDENT}${description[itemindex]}
   echo $FULLLINE
done
unset IFS
}

function usedValues() {
IFS='%'  # this is required because bash automatically strips whitespace otherwise
INDENT=40
echo ""
echo "Values to be used:"
for index in ${!allargedoptions[*]}
do
   itemindex=${allargedoptions[$index]}
   FLAG="    "${options[itemindex]}"                                 "
   FULLLINE=${FLAG:0:40}${argoptionsset[itemindex]}
   echo $FULLLINE
done
unset IFS
}

function getoptionarg(){
   local _findoption=$1
   local _result=""
   shift
   for word in $@
      do
         searchstring="^--""$_findoption""="
         if [[ $word =~ $searchstring ]]
         then
            local _result=${word:((${#searchstring}-1))}
         fi
      done

   echo $_result
}


function findoption(){
   local _findoption=$1
   shift
   for word in $@
      do
      searchstring="^--""$_findoption""="
      if [[ $word =~ $searchstring ]]
      then
         return 0
      fi
      searchstring="--""$_findoption"
      if [[ $word == $searchstring ]]
      then
         return 0
      fi
      done

   return 1
}

function evalNoOptArgs(){
   for option in "${takesnoargument[@]}";
   do
      if findoption ${options[$option]} $@; then
         eval ${options[$option]}Set=1     # set this flag
         noargoptionsset[${#noargoptionsset[*]}]=${options[$option]}
      fi
   done
}

function getMaxVarLength(){
local lenVar=0
for item in $@
do 
   lenItem=${#item}                   # get max name length
   if [ $lenItem -gt $lenVar ]; then
      lenVar=$lenItem
   fi
done
echo $lenVar
}
