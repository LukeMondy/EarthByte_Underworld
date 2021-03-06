import sys, os, subprocess, platform

EnsureSConsVersion(0, 98)

# Colours
colours = {}
colours['cyan']   = '\033[96m'
colours['purple'] = '\033[95m'
colours['blue']   = '\033[94m'
colours['green']  = '\033[92m'
colours['yellow'] = '\033[93m'
colours['red']    = '\033[91m'
colours['end']    = '\033[0m'

# If the output is not a terminal, remove the colours.
if not sys.stdout.isatty():
   for key, value in colours.iteritems():
      colours[key] = ''

compile_source_message = '%sCompiling %s==> %s$SOURCE%s' % \
   (colours['blue'], colours['purple'], colours['yellow'], colours['end'])

link_program_message = '%sLinking Program %s==> %s$TARGET%s' % \
   (colours['cyan'], colours['purple'], colours['yellow'], colours['end'])

link_library_message = '%sLinking Library %s==> %s$TARGET%s' % \
   (colours['cyan'], colours['purple'], colours['yellow'], colours['end'])

ranlib_library_message = '%sRanlib Library %s==> %s$TARGET%s' % \
   (colours['cyan'], colours['purple'], colours['yellow'], colours['end'])

install_message = '%sInstalling %s==> %s$TARGET%s'% \
   (colours['green'], colours['purple'], colours['yellow'], colours['end']) 

# CUSTOMISE THE ENVIRONMENT HERE.
env = Environment(ENV=os.environ,
                  tools=['default', 'pcutest', 'stg', 'dist', 'doc'],
                  toolpath=['StGermain/pcu/script', 'StGermain/script', 'script'])

# Check if scons is launched with detail flag.
detail = ARGUMENTS.get('detail', 0)
Help("""
SCons Build Options:
    Type: './scons.py detail=1' to build Underworld showing full detail non-coloured stdout" 
""" )

# If detail flag is not set, let's use our own build strings.
if not int(detail):
    env['CXXCOMSTR'] = compile_source_message
    env['CCCOMSTR'] = compile_source_message
    env['SHCCCOMSTR'] = compile_source_message
    env['SHCXXCOMSTR'] = compile_source_message
    env['ARCOMSTR'] = link_library_message
    env['RANLIBCOMSTR'] = ranlib_library_message
    env['SHLINKCOMSTR'] = link_library_message
    env['LINKCOMSTR'] = link_program_message
    env['INSTALLSTR'] = install_message

# Load CREDO, the system testing tool
env.Tool('credosystest', toolpath=['credo/scons'])

# Needed for Darwin.
env['_abspath'] = lambda x: File(x).abspath

# Ludicrous-speed!
env.Decider("MD5-timestamp")

# Load configuration.
values = {}
if(not os.path.isfile('config.cfg')):
   print 'Error: \'config.cfg\' file not found.  The configuration script (configure.py) must be run successfully before you can build (via scons.py).'
   sys.exit()

execfile("config.cfg", globals(), values)
env._dict.update(values)

# Set LIBPATH and RPATH for libraries
if env['prefix'] == env.GetLaunchDir():
   env['RPATH'] += [os.path.join(env['prefix'], env['build_dir'], "lib")]
   env['LIBPATH'] = [os.path.join(env['prefix'], env['build_dir'], "lib")] + env['LIBPATH']
else:
   env['RPATH'] += [os.path.join(env['prefix'], "lib")]
   env['LIBPATH'] = [os.path.join(env['prefix'], "lib")] + env['LIBPATH']

# Check if we're using 64bit.
if platform.architecture()[0] == '64bit':
    env.AppendUnique(CPPDEFINES=[('SYSTEM_SIZEOF_LONG', 8)])

# Need to manipulate the build directory to keep SCons happy. Because of SCons' target
# rules we need to make the build directory a default target.
env["build_dir"] = os.path.join(env.GetLaunchDir(), env["build_dir"])
env["prefix"] = os.path.join(env.GetLaunchDir(), env["prefix"])
env["INST_BUILD_DIR"] = env["build_dir"]
env["INST_PREFIX"] = env["prefix"]
env.Default(env["build_dir"])

# Add the build directory's include path.
env.AppendUnique(CPPPATH=env['build_dir'] + '/include')

# Need to define the extension for shared libraries as well
# as the library directory.
ext = env['ESCAPE']('"' + env['SHLIBSUFFIX'][1:] + '"')
lib_dir = env['ESCAPE']('"' + os.path.abspath(env['build_dir']) + '/lib' + '"')
env.AppendUnique(CPPDEFINES=[('MODULE_EXT', ext), ('LIB_DIR', lib_dir)])

# Include the library path.
env.AppendUnique(LIBPATH=env['build_dir'] + '/lib')
env.AppendUnique(RPATH=env.Dir(env['build_dir'] + '/lib').abspath)

# If we have no shared libraries, include a pre-processor definition to
# prevent modules from trying to load dynamically.
if not env['shared_libs']:
    env.AppendUnique(CPPDEFINES=['NOSHARED'])

# Need to insert some 'HAVE_*' definitions based on what packages we
# found during configuration.
if 'HAVE_HDF5' in env['CPPDEFINES']:
    env.AppendUnique(CPPDEFINES=["READ_HDF5", "WRITE_HDF5"])

# If we were given a prefix other than the default, tell StGermain where to
# find XML include files.
if env['prefix'] != env.GetLaunchDir():
    env.AppendUnique(CPPDEFINES=[('STG_INCLUDE_PATH', env['ESCAPE']('"' + env['prefix'] + '/lib/StGermain"'))])

# Make sure 'install' has a proper target.
env.Alias("install", env["prefix"])

# INSERT TARGETS HERE.

# Make a copy of the config.log and config.cfg in the build_dir
dir = Dir(env['build_dir']).abspath
if not os.path.exists(dir):
   os.makedirs(dir)

env.Install(env['build_dir'], 'config.log')
env.Install(env['build_dir'], 'config.cfg')
env.Install(os.path.join(env['build_dir'],'share/StGermain/scripts'), Glob('script/*.py*'))
if env['prefix'] != env.GetLaunchDir():
   env.Install(env['prefix'], 'config.log')
   env.Install(env['prefix'], 'config.cfg')
   env.Install(os.path.join(env['prefix'],'share/StGermain/scripts'), Glob('script/*.py*'))

Export('env')


#########################
# Setup hg precommit hook
# to strip .ipynb
#########################
hook_found=False
pos=-1
lineCount=0
# first look if hook already exists
try:
   hgrc = open("../.hg/hgrc", "r+")
except IOError:
   pass
else:
   line=hgrc.readline()
   while line:
     if line.find("pre-commit=\"./libUnderworld/script/precommit.py\"",0) == 0:
       print "Found precommit hook declaration"
       hook_found=True

     # if hooks specified find the position in file
     if line.find("[hooks]",0) == 0:
       pos = lineCount 
       
     lineCount=lineCount+1
     line=hgrc.readline() # iterate

   hgrc.close()

   adding=""
# if hook not found, add the magic in the .hg/hgrc
   if hook_found == False:
     print "Adding magic precommit hook delaration in .hg/hgrc"

     # read in contents of hgrc
     hgrc = open("../.hg/hgrc", "r")
     contents=hgrc.readlines()
     hgrc.close()

     # find place to add and create directive
     if pos == -1:
       pos=len(contents)
       adding = "\n[hooks]\n"

     adding += "pre-commit=\"./libUnderworld/script/precommit.py\"\n"
     # insert directive
     contents.insert(pos+1, adding)

     # replace new contents with old
     hgrc=open("../.hg/hgrc","w")
     contents = "".join(contents)
     hgrc.write(contents)
     hgrc.close()

############################
# Finished hg precommit hook
############################


SConscript('StGermain/SConscript',
           variant_dir=env['build_dir'] + '/StGermain',
           duplicate=0)
env.Prepend(LIBS=['pcu'])
env.Prepend(LIBS=['StGermain'])

SConscript('StgDomain/SConscript',
           variant_dir=env['build_dir'] + '/StgDomain',
           duplicate=0)
env.Prepend(LIBS=['StgDomain'])

SConscript('StgFEM/SConscript',
           variant_dir=env['build_dir'] + '/StgFEM',
           duplicate=0)
env.Prepend(LIBS=['StgFEM'])

SConscript('PICellerator/SConscript',
           variant_dir=env['build_dir'] + '/PICellerator',
           duplicate=0)
env.Prepend(LIBS=['PICellerator'])

SConscript('Underworld/SConscript',
           variant_dir=env['build_dir'] + '/Underworld',
           duplicate=0)
env.Prepend(LIBS=['Underworld'])

if env['with_importers']:
    SConscript('ImportersToolbox/SConscript',
               variant_dir=env['build_dir'] + '/ImportersToolbox',
               duplicate=0)
    env.Prepend(LIBS=['ImportersToolbox'])

if env['with_solvers']:
    SConscript('Solvers/SConscript',
               variant_dir=env['build_dir'] + '/Solvers',
               duplicate=0)
    env.Prepend(LIBS=['Solvers'])

if env['with_spherical']:
    SConscript('Spherical/SConscript',
               variant_dir=env['build_dir'] + '/Spherical',
               duplicate=0)
    env.Prepend(LIBS=['Spherical'])

if env['with_viscoelastic']:
    SConscript('Viscoelastic/SConscript',
               variant_dir=env['build_dir'] + '/Viscoelastic',
               duplicate=0)
    env.Prepend(LIBS=['Viscoelastic'])




# Dump package config.
filename = env['build_dir'] + '/lib/pkgconfig/stgermain.pc'

env.Dist("underworld-%s"%env.GetOption("dist_version"),
         ["configure.py", "SConstruct", "config", "script", "StGermain",
          "StgDomain", "StgFEM", "PICellerator", "Underworld",])

# NB: help() printout about testing that used to be here moved to credo/scons.

    
#Run the functions to create the Alias commands for the documentation of stgUnderworld. 


env.Alias("doc-codex", None, env.Action(env.AddCodexSuite()))
env.Alias("doc", None, env.Action(env.AddAllSuite()))
env.Alias("doc-doxygen", None, env.Action(env.AddDoxygenSuite()))
env.Alias("doc-doxygenlite", None, env.Action(env.AddDoxygenLiteSuite()))

# TODO: test targets and master suite runners currently go below,
# see CREDO tool.
# Ideally if could do all target stuff properly, this should go back in the
#  CREDO tool.
Import('LOWRES_SUITES')
Import('INTEGRATION_SUITES')
Import('VISUALISATION_SUITES')
Import('CONVERGENCE_SUITES')
Import('SCIBENCH_SUITES')
lowresSuiteRun = env.RunSuites( 
    Dir(os.path.join(env['TEST_OUTPUT_PATH'], env["CHECK_LOWRES_TARGET"])),
    LOWRES_SUITES)
env.AlwaysBuild(lowresSuiteRun)
env.Alias(env["CHECK_LOWRES_TARGET"], lowresSuiteRun)
intSuiteRun = env.RunSuites(
    Dir(os.path.join(env['TEST_OUTPUT_PATH'], env["CHECK_INTEGRATION_TARGET"])),
    INTEGRATION_SUITES)
env.AlwaysBuild(intSuiteRun)
env.Alias(env["CHECK_INTEGRATION_TARGET"], intSuiteRun)
vizSuiteRun = env.RunSuites(
    Dir(os.path.join(env['TEST_OUTPUT_PATH'], env["CHECK_VISUALISATION_TARGET"])),
    VISUALISATION_SUITES)
env.AlwaysBuild(vizSuiteRun)
env.Alias(env["CHECK_VISUALISATION_TARGET"], vizSuiteRun)
cvgSuiteRun = env.RunSuites(
    Dir(os.path.join(env['TEST_OUTPUT_PATH'], env["CHECK_CONVERGENCE_TARGET"])),
    CONVERGENCE_SUITES)
env.AlwaysBuild(cvgSuiteRun)
env.Alias(env["CHECK_CONVERGENCE_TARGET"], cvgSuiteRun)
scibenchSuiteRun = env.RunSuites( 
    Dir(os.path.join(env['TEST_OUTPUT_PATH'], env["CHECK_SCIBENCH_TARGET"])),
    SCIBENCH_SUITES)
env.AlwaysBuild(scibenchSuiteRun)
env.Alias(env["CHECK_SCIBENCH_TARGET"], scibenchSuiteRun)
# Run the lowres checks as part of default and complete
env.Alias(env['CHECK_DEFAULT_TARGET'], env['CHECK_LOWRES_TARGET'])
env.Alias(env['CHECK_COMPLETE_TARGET'], env['CHECK_LOWRES_TARGET'])
# For the others, just add to the complete target
env.Alias(env['CHECK_COMPLETE_TARGET'], env['CHECK_INTEGRATION_TARGET'])
env.Alias(env['CHECK_COMPLETE_TARGET'], env['CHECK_CONVERGENCE_TARGET'])
#Don't run scibench as part of 'complete', since could run for 
# -very- long time.
#env.Alias(env['CHECK_COMPLETE_TARGET'], env['CHECK_SCIBENCH_TARGET'])

