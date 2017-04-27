import numpy as np
#from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np

filelist=[]

for i in range(200,805,10):
    filelist.append("Stats_%s.txt" %i)

fig, ax = plt.subplots()
ax.grid(True)
#ax.axis([0.,500.,100.,850.])
ax.axis([0.,500.,0.,1.5])
initialTemp = 200.
for fname in filelist:
#    print fname
    data=np.loadtxt(fname)
    x=data[:,0]
    y=data[:,2]
#    ax.plot(x,y)
    ax.plot(x,y/initialTemp)
    initialTemp += 10

plt.xlabel('Step Number')
plt.ylabel('Temperature [K]')
plt.show()
