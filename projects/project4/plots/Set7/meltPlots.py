import numpy as np
#from cycler import cycler
import matplotlib.pyplot as plt

filelist=[]

for i in range(560,590,5):
    filelist.append("Stats_%s.txt" %i)

print "File list generated"

#fig, (ax0,ax1) = plt.subplots(nrows=2)
fig, (ax0) = plt.subplots()
ax0.grid(True)
#ax1.grid(True)
ax0.axis([0.,10000.,0.,0.005])
#ax1.axis([200.,1000.,0.,0.004])
initialTemp = 560.
for fname in filelist:
    eqTemp = initialTemp/2.
    print fname
    data=np.loadtxt(fname)
    x=data[:,0]
    y=data[:,6]
    lname = 'T = %s' %eqTemp
    ax0.plot(x,y, label=lname)
#    ax1.plot(x,y2)
    initialTemp += 5

plt.xlabel(r'Step Number ($i$)')
#ax0.set_ylabel('Temperature [K]')
ax0.set_ylabel(r'$D$')
plt.legend(loc=1)
plt.show()
