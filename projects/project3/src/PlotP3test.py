import os,sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileEuler  = "../Benchmark/SEsystemEuler.out"
fileVerlet = "../Benchmark/SEsystemVerlet.out"

data = np.loadtxt(fileEuler)
x = data[0,:]
y = data[1,:]

dataV = np.loadtxt(fileVerlet)
xV = dataV[0,:]
yV = dataV[1,:]

plotE = plt.plot(x,y,'b')
plotV = plt.plot(xV,yV,'r')
plt.axis([-2.,2.,-2.,2.])
plt.show()


