#Fock diagram

from scipy import *
import matplotlib.pyplot as plt
import numpy as np

def Sigma_x(EF, k, lam):
    kF = sqrt(EF)
    pp = 0
    if lam>0:
        pp = lam/kF*(arctan((k+kF)/lam)-arctan((k-kF)/lam))
    qq = 1 - pp - (lam**2+kF**2-k**2)/(4*k*kF)*log((lam**2+(k-kF)**2)/(lam**2+(k+kF)**2))
    return -2*kF/(pi)*qq

EF=1.91
lam=0.1
k=np.linspace(0.001, 20, 1000)
Sigma=Sigma_x(EF, k, lam)
shift=Sigma_x(EF, EF, lam)
Sigma-=shift
for i in range(len(k)):
    print k[i], Sigma[i]
plt.plot(k, Sigma)
plt.show()
