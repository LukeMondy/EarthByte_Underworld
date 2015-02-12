# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from libUnderworld import *

# <codecell>

from libUnderworld import c_arrays

# <codecell>

location = c_arrays.DoubleArray(2)
velocity = c_arrays.DoubleArray(2)
pressure = c_arrays.DoubleArray(1)
stress = c_arrays.DoubleArray(3)
strain_rate = c_arrays.DoubleArray(3)

# <codecell>

location[0]=0.7
location[1]=0.7

# <codecell>

analytic.solcx(location.cast(), 10000.0,1.0, 0.5, 2, velocity.cast(), pressure.cast(), None, None) # <- works fine with None's
#solcx(location, eta_A, eta_B, x_c, n, velocity, pressure, stress, strain_rate)
analytic.solcx(location.cast(), 1000000.0,1.0, 0.55, 2, velocity.cast(), pressure.cast(), stress.cast(), strain_rate.cast())

# <codecell>

print velocity[0], velocity[1], pressure[0], stress[0], stress[1], stress[2], strain_rate[0], strain_rate[1], strain_rate[2]

# <codecell>

#solkx(location, sigma, m, n, B, velocity, pressure, stress, strain_rate, None)
analytic.solkx(location.cast(), 1.0, 2, 3, 1.686, velocity.cast(), pressure.cast(), stress.cast(), strain_rate.cast(), None)

# <codecell>

print velocity[0], velocity[1], pressure[0], stress[0], stress[1], stress[2], strain_rate[0], strain_rate[1], strain_rate[2]

# <codecell>

#solkz(location, sigma, m, n, B, velocity, pressure, stress, strain_rate)
analytic.solkz(location.cast(), 1.0, 2, 3, 1.686, velocity.cast(), pressure.cast(), stress.cast(), strain_rate.cast())

# <codecell>

print velocity[0], velocity[1], pressure[0], stress[0], stress[1], stress[2], strain_rate[0], strain_rate[1], strain_rate[2]

# <codecell>


