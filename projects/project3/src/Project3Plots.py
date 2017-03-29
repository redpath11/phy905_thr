import os,sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

filein     = "../Benchmark/SEsystemEuler.out"
fileV      = "../Benchmark/SEsystemVerlet.out"
file3BE    = "../Benchmark/SJEsystemEuler.out"
fileVec3   = "../Benchmark/vec3test.out"
filePlanet = "../Benchmark/PlanetVerlet.out"
fileSun    = "../Benchmark/SunVerlet.out"
fileEarth  = "../Benchmark/EarthVerlet.out"
fileSS     = "../Benchmark/ooSystem.out"

#data = np.loadtxt(filein)
#x = data[0,:]
#y = data[1,:]

#dataV = np.loadtxt(fileV)
#xV = dataV[0,:]
#yV = dataV[1,:]

#data = np.loadtxt(file3BE)
#xE = data[0,:]
#yE = data[1,:]
#xJ = data[4,:]
#yJ = data[5,:]

#data = np.loadtxt(fileSun)
#xS = data[:,0]
#yS = data[:,1]
#zS = data[:,2]

#data = np.loadtxt(fileEarth)
#xE = data[:,0]
#yE = data[:,1]
#zE = data[:,2]

data = np.loadtxt(fileSS)

xV = data[:,0]
yV = data[:,1]
zV = data[:,2]

xE = data[:,3]
yE = data[:,4]
zE = data[:,5]

xM = data[:,6]
yM = data[:,7]
zM = data[:,8]

xJ = data[:,9]
yJ = data[:,10]
zJ = data[:,11]

xSt = data[:,12]
ySt = data[:,13]
zSt = data[:,14]

xS = data[:,15]
yS = data[:,16]
zS = data[:,17]
#plotE = plt.plot(x,y,'b')
#plotV = plt.plot(xV,yV,'r')
#plotE = plt.plot(xE,yE,'b');
#plotJ = plt.plot(xJ,yJ,'r');
#plotE = plt.plot(x,y,'r')
#plt.show()
fig = plt.figure()
ax1 = fig.add_subplot(111,projection='3d')
ax1.plot(xS,yS,zS,'y')
ax1.plot(xV,yV,zV,'y')
ax1.plot(xE,yE,zE,'b')
ax1.plot(xM,yM,zM,'r')
ax1.plot(xJ,yJ,zJ,'r')
ax1.plot(xSt,ySt,zSt,'m')
plt.show()
