#!/usr/bin/env python

'''
Simple 1d diffusion using crank-nicolson
'''

import scipy.linalg
import numpy
import time
import matplotlib.pyplot as plt

UWcouple = True
nx_step=80
diffusivity = 0.00025
dx = 1./float(nx_step)
dt = 0.1*dx**2/diffusivity
t = 0.
sync_interval = 5.
# lets make a nicer timestep size
steps_per_sync = int(sync_interval/dt)
steps_per_sync+=1
dt = sync_interval/float(steps_per_sync)
r =  0.5*diffusivity*dt/(dx)**2
maxt = 20.


coefa = numpy.zeros( (nx_step, nx_step) )
u     = numpy.zeros( nx_step )
uwin  = numpy.zeros( nx_step )
b     = numpy.zeros( nx_step )

# set init condition
u[:]=0.4
for j in numpy.arange(0+30, nx_step-30):
    u[j]=0.6

u0 = u.copy()

print 'Running 1d diffusion'
if UWcouple:
    # lets write all to a file
    f = open('SPOuput.txt', 'w')
    for j in numpy.arange(0, nx_step):
        outstr = str(j) +' '+ str(u[j]) + '\n'
        f.write(outstr)
    f.close()
    print 'outputted SP surface to SPOutput.txt'

    # lets open the maestro file, and set it to U
    f = open('maestro', 'w')
    f.write('U')
    f.close()
    print 'Waiting for UW to solve until sync time of', sync_interval

for j in numpy.arange(0, nx_step):
    coefa[j, (j-1)%nx_step] = - r
    coefa[j, (j  )%nx_step] = 1. + 2.*r
    coefa[j, (j+1)%nx_step] = - r

inva = scipy.linalg.inv(coefa)         # get inverse of A

nextsync = t + sync_interval
wait = False
if UWcouple:
    wait = True

#plt.plot(u)
#plt.plot(u0)
#plt.ylim([0,1])
#plt.show()

while t < maxt:
    #---- main loop ----
    # lets open the maestro file, and wait for our turn
    while wait:
        f = open('maestro', 'r')
        if f.read() == 'L':
            wait = False
            print 'Ok, UW at sync time.  Now run SP.'
            # ok, lets read in uw file
            data = numpy.genfromtxt("uw_output.ascii")
            for datum in data:
                uwin[int(datum[0])] = float(datum[1])
        else:
            time.sleep(2)
        f.close()
        
    
    for j in numpy.arange(0, nx_step):
        b[j] =    (      r)*u[(j+1)%nx_step]  \
                + (1.-2.*r)*u[(j  )%nx_step]  \
                + (      r)*u[(j-1)%nx_step]

    u = numpy.dot(inva,b)
    # ok, lets slowly introduce UW deformations
    u[:] += uwin[:]/float(steps_per_sync)
    t = t + dt
    #plt.plot(u)
    #plt.plot(u0)
    #plt.ylim([0,1])
    #plt.show()

    print 'time =', t
    if t>=nextsync and UWcouple:
        #plt.plot(u)
        #plt.plot(u0)
        #plt.ylim([0,1])
        #plt.show()
        
        wait = True
        print 'SP at next synctime, passing back to UW.  Current time is: ', t
        nextsync = t + sync_interval
        # lets write all to a file
        f = open('SPOuput.txt', 'w')
        for j in numpy.arange(0, nx_step):
            outstr = str(j) +' '+ str(u[j]) + '\n'
            f.write(outstr)
        f.close()
        # update maestro
        f = open('maestro', 'w')
        f.write('U')
        f.close()

 
