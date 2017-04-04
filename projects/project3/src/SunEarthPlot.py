import os,sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileEuler  = "../Benchmark/ooEuler.out"
fileVerlet = "../Benchmark/ooSystem.out"

data = np.loadtxt(fileEuler)
xE = data[:,0]
yE = data[:,1]
zE = data[:,2]

xS = data[:,3]
yS = data[:,4]
zS = data[:,5]

data = np.loadtxt(fileVerlet)
xEv = data[:,0]
yEv = data[:,1]
zEv = data[:,2]

xSv = data[:,3]
ySv = data[:,4]
zSv = data[:,5]

#fig = plt.figure()
#ax1 = fig.add_subplot(111,projection='3d')
#ax1.plot(xE,yE,zE,'b')
#ax1.plot(xS,yS,zS,'y')
#plt.show()

fig, ax = plt.subplots()
plotE = ax.plot(xE,yE,'b',label='Earth (Euler)')
plotS = ax.plot(xS,yS,'y*',markersize=10.,label='Sun')
plotEv = ax.plot(xEv,yEv,'r',label='Earth (VVerlet)')
ax.grid(True)
#plotSv = plt.plot(xSv,ySv,'y*',markersize=10.)
ax.axis([-50.,50.,-50.,50.])
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title('Sun-Earth System')
plt.legend(loc=1)
plt.show()
