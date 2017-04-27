import numpy as np
import matplotlib.pyplot as plt

filename = "../src/statistics.txt"

data=np.loadtxt(filename)
# plot temperature vs. time step
x=data[:,0]
y=data[:,2]
initialTemp = 300.

fig, ax = plt.subplots()
ax.grid(True)
ax.axis([0,1000.,0,1.5])
ax.plot(x,y/initialTemp)
plt.show()
