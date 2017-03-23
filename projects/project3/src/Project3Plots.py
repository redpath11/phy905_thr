import os,sys
import numpy as np
from matplotlib import pyplot as plt

filein = "../Benchmark/SEsystemEuler.out"
fileV  = "../Benchmark/SEsystemVerlet.out"
file3BE= "../Benchmark/SJEsystemEuler.out"
fileVec3 = "../Benchmark/vec3test.out"

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

data = np.loadtxt(fileVec3)
x = data[:,0]
y = data[:,1]

#plotE = plt.plot(x,y,'b')
#plotV = plt.plot(xV,yV,'r')
#plotE = plt.plot(xE,yE,'b');
#plotJ = plt.plot(xJ,yJ,'r');
plotE = plt.plot(x,y,'r')
plt.show()
