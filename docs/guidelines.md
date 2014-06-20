

Design
======
User Interface?
Naming convections - for classes, instances, input parameters?
Class overview?


Directory Structure
===================

* libUnderworld            - contains backend stack of tools (Underworld, StgFEM, etc)
 * libUnderworldPy         - contains swig generated wrappers for backend as well as auxiliary implementations (such as StGermain_Tools)
* underworld               - python only directory structure containing python frontend routines and supporting scripts
* docs                     - contains various documentation including this file
* InputFiles               - contains various examples input files



Coding Style
============

Whitespace
----------
PEP 8
http://legacy.python.org/dev/peps/pep-0008/

Most importantly:
Spaces should be used for tabs.
Indent levels should be 4 spaces deep.

autopep8 tool can be used to tidy when necessary:
https://pypi.python.org/pypi/autopep8/


* Comments?
* Error handling?
 * Reporting?
 * Checking / exceptions?


License 
=======



Development
===========
* Versions?
* Development Workflow?
* Bug fixes & reporting?
* New Features / Redesigned interface or API?
* Release candidates?
* Distribution?

Testing
=======
* Unit tests?
* Regression tests?


Useful Link
===========
http://matplotlib.org/devel/gitwash/git_development.html
https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
http://docs.scipy.org/doc/numpy/dev/index.html
