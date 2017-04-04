import os,sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileEuler  = "../Benchmark/ooEuler.out"
fileVerlet = "../Benchmark/ooSystem.out"

#data = np.loadtxt(fileEuler)
data = np.loadtxt(fileVerlet)
xM = data[:,0]
yM = data[:,1]
zM = data[:,2]

xS = data[:,3]
yS = data[:,4]
zS = data[:,5]

#fig = plt.figure()
#ax1 = fig.add_subplot(111,projection='3d')
#ax1.plot(xM,yM,zM,'b')
#ax1.plot(xS,yS,zS,'y')
plotM = plt.plot(xM,yM,'m',label='Mercury')
plotS = plt.plot(xS,yS,'y*',markersize=10.,label='Sun')
plt.axis([-1.,1.,-1.,1.])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Mercury Perihelion Precession')
plt.legend(loc=1)

plt.show()

#plotE = plt.plot(xE,yE,'b');
#plotS = plt.plot(xJ,yJ,'r');
