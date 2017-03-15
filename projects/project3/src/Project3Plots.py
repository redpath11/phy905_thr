import os,sys
import numpy as np
from matplotlib import pyplot as plt

filein = "../Benchmark/SEsystemEuler.out"
fileV  = "../Benchmark/SEsystemVerlet.out"

data = np.loadtxt(filein)
x = data[0,:]
y = data[1,:]

dataV = np.loadtxt(fileV)
xV = dataV[0,:]
yV = dataV[1,:]

plotE = plt.plot(x,y,'b')
plotV = plt.plot(xV,yV,'r')
plt.show()
