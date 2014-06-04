#!/usr/bin/env python
import sys, subprocess, shutil, os

dir = "../libUnderworld/build"

# pull the config directory from the StGermain build
subp = subprocess.Popen(
	'hg clone https://www.underworldproject.org/hg/SConfigure config ', shell=True
)
subp.wait()

sconsBin = os.path.join('config', 'scons', 'scons.py')

subp = subprocess.Popen(
	sconsBin + ' --config=force -f SConfigure ' + 
	'--build-dir=build ' +
	'--pcu-dir=${UW_DIR} ' +
	'--stgermain-dir=${UW_DIR} ' +
	'--stgdomain-dir=${UW_DIR} ' +
	'--stgfem-dir=${UW_DIR} ' +
	'--picellerator-dir=${UW_DIR} ' +
	'--glucifer-dir=${UW_DIR} ' +
	'--underworld-dir=${UW_DIR} ' + ' '.join(sys.argv[1:]), shell=True
)
subp.wait()
