#!/usr/local/macports/bin/python
from numpy import loadtxt
import pylab as plt

data = loadtxt('build/derp.out', delimiter=' ', dtype=float)
print data

rad      = data[:,  0]
theta1   = data[:,  1]
theta2   = data[:,  2]
theta3   = data[:,  3]
theta127 = data[:,127]

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)

ax.plot(theta1  , rad)
ax.plot(theta2  , rad)
ax.plot(theta3  , rad)
ax.plot(theta127, rad)

plt.show()
