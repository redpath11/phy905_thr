import os,sys
import numpy as np
from matplotlib import pyplot as plt

filein = "../Benchmark/SEsystemEuler.out"

data = np.loadtxt(filein)
x = data[0,:]
y = data[1,:]

plot = plt.plot(x,y)
plt.show()
