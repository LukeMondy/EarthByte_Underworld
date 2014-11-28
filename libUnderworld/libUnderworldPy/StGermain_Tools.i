/* -*- C -*-  (not really, but good for syntax highlighting) */

%module StGermain_Tools

%import "StGermain.i"

%{
/* Includes the header in the wrapper code */
#include "StGermain_Tools.h"
#include <string.h>
#include <StGermain/StGermain.h>

%}

/* Parse the header file to generate wrappers */
%include <argcargv.i>
%apply (int ARGC, char **ARGV) { (int argc, char *argv[]) }
%include "StGermain_Tools.h"
