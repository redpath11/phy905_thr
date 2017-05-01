import numpy as np
import matplotlib.pyplot as plt

filename = "../src/statistics.txt"
#filename = "Set3/Stats_50.txt"

data=np.loadtxt(filename)
# plot temperature vs. time step
#x=data[:,0]
#y=data[:,2]
#initialTemp = 300.

# plot diffusion const vs. time step
x=data[:,0]
y=data[:,6]

# plot diffusion const vs. temp
#x=data[:,2]
#y=data[:,6]


fig, ax = plt.subplots()
ax.grid(True)
#ax.axis([0,1000.,0,1.5])
ax.axis([0,10000.,0,0.025])
#ax.axis([0,600.,0,0.05])
ax.plot(x,y)
#ax.plot(x,y/initialTemp)
plt.show()
