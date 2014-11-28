import os, sys, subprocess

from runipy.notebook_runner import NotebookRunner, NotebookError
from IPython.nbformat.current import read

# create the test directory if needed
dir = "./testDir/"
try:
   os.stat(dir)
except:
   os.mkdir(dir)

command = "runipy -h"
try:
  subprocess.check_call(command.split())
except:
  print "\nCould execute test because I can't execute 'runipy'"
  print "Make sure 'runipy' is installed"
  print "$ pip install runipy\n\n"
  sys.exit(1)


command = "runipy -q"

nfails=0
list_fails=[]
logFile = open(dir+"testing.log", "w") # create test log file

for f in os.listdir('.'):
  if f.endswith(".ipynb"):

    fname = f             # file name

    exe = command.split() # make command a list
    exe.append(f)         # append filename

    outfile = dir + f
    exe.append(outfile)   # append output

    logFile.write("EXECUTION of " + fname); logFile.flush()

# try run runipy on given notebook
    try:
      subprocess.check_call( exe )
    except subprocess.CalledProcessError:
      logFile.write(" .... ERROR (see "+outfile+")\n\n"); logFile.flush()
      nfails = nfails+1
      list_fails.append(fname)
    else:
      logFile.write(" .... PASS \n\n"); logFile.flush()
       

# Report in testing.log
logFile.write("Number of fails " + str(nfails) +":\n"); logFile.flush()
logFile.write(str(list_fails)); logFile.flush()

# Report to stdout
print "\n\n Total: Number of fails " + str(nfails)  
print list_fails

logFile.close()

if nfails > 0:
   sys.exit( nfails )
else:
   sys.exit(0)
