import os,sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileVerlet = "../Benchmark/ooSEJ.out"

data = np.loadtxt(fileVerlet)
xE = data[:,0]
yE = data[:,1]
zE = data[:,2]

xJ = data[:,3]
yJ = data[:,4]
zJ = data[:,5]

xS = data[:,6]
yS = data[:,7]
zS = data[:,8]

#fig = plt.figure()
#ax1 = fig.add_subplot(111,projection='3d')
#ax1.plot(xE,yE,zE,'b')
#ax1.plot(xJ,yJ,zJ,'r')
#ax1.plot(xS,yS,zS,'y')
plotS = plt.plot(xS,yS,'y*',markersize=10.,label='Sun')
plotE = plt.plot(xE,yE,'b',label='Earth')
plotJ = plt.plot(xJ,yJ,'r',label='Jupiter')
plt.axis([-6.,6.,-6.,6.])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Planet Orbits - Three body system')
plt.legend(loc=1)
plt.show()

#plotE = plt.plot(xE,yE,'b');
#plotS = plt.plot(xJ,yJ,'r');
