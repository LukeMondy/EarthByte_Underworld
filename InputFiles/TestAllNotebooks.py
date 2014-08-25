import os, sys, subprocess

# create the test directory if needed
dir = "./testDir/"
try:
   os.stat(dir)
except:
   os.mkdir(dir)

# turn all ipython notebooks into python scripts under ./testDir/
command = "ipython nbconvert --to python --stdout "

for f in os.listdir('.'):
  if f.endswith(".ipynb"):
      outfile = dir + f[:-6] + ".py"
      exe = command.split() # make command a list
      exe.append(f)         # append filename
      print exe
      with open(outfile, "w") as outfile:
         retCode = subprocess.call(exe, stdout=outfile)
      if retCode != 0:
         print "WTF nbconvert can't convert " + f + " ... exiting"
         sys.exit(1)

# change directory to ./testDir
os.chdir(dir)

nfails=0
list_fails=[]
logFile = open("testing.log", "w")

# run each .py and record it error code. If nonzero assume failure
for f in os.listdir('.'):
   if f.endswith(".py"):
      logFile.write("\n\n\n ***********\nEXECUTION of " + f+"\n"); logFile.flush()

      retCode = subprocess.call(['python', f], stdout=logFile, stderr=logFile )

      if retCode != 0:
        logFile.write(" **** ERROR ****\n"); logFile.flush()

        nfails = nfails + 1
        list_fails.append(f)
   

print "Number of fails " + str(nfails)  
print list_fails

if nfails > 0:
   sys.exit( nfails )
else:
   sys.exit(0)
