#!/usr/local/macports/bin/python
from numpy import loadtxt
import pylab as plt

mu, s, rad = loadtxt('build/derp.out', unpack=True)

fig = plt.figure()
#ax = fig.add_subplot(111)
ax = fig.add_subplot(111, polar=True)

ax.plot(mu, rad, 'x')
#ax.set_rmin(6.0e14)
#ax.set_rmax(1.0e15)

plt.show()
