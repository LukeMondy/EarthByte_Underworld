

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

This command seems to give good results (not overly aggressive):
autopep8 -v -i -r  --ignore E201,E202,E501,E221,E251 .

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

Jenkins, CI server
------------------
We run a contiguous integration system Jenkins. Currently it downloads the code each night, configures, compiles and runs the units and system tests. If an error is detected in any of these processes the code is labelled as a Failure.

See https://130.56.248.95:8080/job/underworld2/

To test ipython notebooks perhaps we can use scripts like this https://gist.github.com/shoyer/7497853


Useful Link
===========
http://matplotlib.org/devel/gitwash/git_development.html
https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
http://docs.scipy.org/doc/numpy/dev/index.html
