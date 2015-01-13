#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import underworld as uw

res=16
n=4

dx  = np.zeros(n)
err = np.zeros(n)

for ii in np.arange(n):
	uw.Init("testSemiLagrangian.xml")
	gdict=uw.dictionary.GetDictionary()
	res = res*2
	gdict['elementResI']=res
	gdict['elementResJ']=res
	uw.dictionary.SetDictionary(gdict)
	uw.Construct()
        uw.RunMainLoop()

	dx[ii]  = 2.0/res
        err[ii] = np.loadtxt( 'error.txt' )

log_x = np.log(dx)
log_e = np.log(err)

y   = np.vstack([log_x,np.ones(len(log_x))]).T
m,c = np.linalg.lstsq(y,log_e)[0]
print 'convergence rate: ',m

m_old = 3.0353768958
if m < m_old - 0.1*m_old:
	print 'ERROR: convergence rate has degraded from ', m_old

plt.loglog(dx,err,'o-')
plt.show()
