import sys, os
from matplotlib import pyplot as plt
import numpy as np

def lj(r,epsilon=1.0, sigma=1.0):
    if r>0.:
        return 4.*epsilon*(sigma**12/r**12 - sigma**6/r**6)
    else:
        return None

def flj(r,epsilon=1.0, sigma=1.0):
    if r>0.:
        return (24.*epsilon/r)*(2.*sigma**12/r**12 - sigma**6/r**6)
    else:
        return None

x = np.arange(0.5,3.0,0.001)
vlj  = np.vectorize(lj)
vflj = np.vectorize(flj)

fig, ax = plt.subplots()
ljplot  = ax.plot(x,vlj(x),'b',label=r'$U(r)$')
fljplot = ax.plot(x,vflj(x),'r',label=r'$F(r)$')
ax.axis([0,3.0,-2.5,0.5])
ax.grid(True)
plt.legend(loc=4)
plt.xlabel(r'$r$')
#plt.ylabel(r'$U(r)$')
plt.show()
