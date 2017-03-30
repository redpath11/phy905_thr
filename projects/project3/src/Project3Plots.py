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

xHg = data[:,0]
yHg = data[:,1]
zHg = data[:,2]

xV = data[:,3]
yV = data[:,4]
zV = data[:,5]

xE = data[:,6]
yE = data[:,7]
zE = data[:,8]

xM = data[:,9]
yM = data[:,10]
zM = data[:,11]

xJ = data[:,12]
yJ = data[:,13]
zJ = data[:,14]

xSt = data[:,15]
ySt = data[:,16]
zSt = data[:,17]

xU = data[:,18]
yU = data[:,19]
zU = data[:,20]

xN = data[:,21]
yN = data[:,22]
zN = data[:,23]

xP = data[:,24]
yP = data[:,25]
zP = data[:,26]

xS = data[:,27]
yS = data[:,28]
zS = data[:,29]
#plotE = plt.plot(x,y,'b')
#plotV = plt.plot(xV,yV,'r')
#plotE = plt.plot(xE,yE,'b');
#plotJ = plt.plot(xJ,yJ,'r');
#plotE = plt.plot(x,y,'r')
#plt.show()
fig = plt.figure()
ax1 = fig.add_subplot(111,projection='3d')
ax1.plot(xS,yS,zS,'y')
ax1.plot(xHg,yHg,zHg,'k')
ax1.plot(xV,yV,zV,'y')
ax1.plot(xE,yE,zE,'b')
ax1.plot(xM,yM,zM,'r')
ax1.plot(xJ,yJ,zJ,'r')
ax1.plot(xSt,ySt,zSt,'m')
ax1.plot(xU,yU,zU,'b')
ax1.plot(xN,yN,zN,'b')
ax1.plot(xP,yP,zP,'k')
plt.show()
