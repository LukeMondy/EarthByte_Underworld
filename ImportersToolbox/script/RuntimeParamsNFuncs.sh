#!/bin/bash
#
# Underworld Importer Geothermal model builder
#
# Copyright (C) 2012, Monash University
#
# Contact: John Mansour

# Define input variables variables. Add new items as required.
# ------------------------------------------------------------

flag="run"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Run model once created"
takesnoargument[${#takesnoargument[*]}]=$item
unset flag; unset item; unset testfntext;

flag="mpirun"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Name of the mpirun command on your system"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]="mpiexec"
unset flag; unset item; unset testfntext;

flag="nprocs"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Number of processors to run job with"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]=1
unset flag; unset item; unset testfntext;

flag="bin"
options[${#options[*]}]=$flag
item=$((${#options[*]}-1))
eval ${options[$item]}=$item
description[$item]="Underworld binary executable, including absolute path"
takesargument[${#takesargument[*]}]=$item
defaultvalue[$item]="${DIR}/../../build/bin/Underworld"
unset flag; unset item; unset testfntext;

