import MCMC as mcmc
import pickle
import matplotlib.pyplot as plt
import numpy as np
from pylab import hist, show

n=4
mean=[1,400,30,50]
sigma=[3,20,5,5]


def fluxFunction(x):
    return 1e3*x**-2.7

cov=np.diag([0.5 for i in range(n)])
for i in range(n):
    if i+1 < n:
        cov[i,i+1]=0.25
        if i-1 >= 0:
            cov[i,i-1]=0.25



realFlux=[fluxFunction(i+1) for i in range(n)]
expectedAfterSmearing=cov.dot(np.array(realFlux))
data_observed=np.random.poisson(expectedAfterSmearing)            

def logp(value):
    theExpectedFlux=cov.dot(np.array(value))
    log=0
    for i in range(n):
        if value[i]<0:
            return -np.inf
        var=data_observed[i]*np.log(theExpectedFlux[i]) - theExpectedFlux[i]
        log+=var

    return log

def generateMCMC(filename):
    a=mcmc.MCMC()
    a.setLogLikelihoodFunction(logp)
    a.
    a.setSteps(10000)
    a.loop()
    f=open(filename,'wb')
    pickle.dump(a,f)
    f.close()
    return a

def loadMCMC(filename):
    f=open(filename,'rb')
    return pickle.load(f)

filename='test.pkl'
a=generateMCMC(filename)
#a=loadMCMC('test.pkl')

burnInLength=600
correlationLength=1
for i in range(n):
    plt.figure()
    plt.subplot(211)

    ar=np.array(a.trace)[:,i][burnInLength::correlationLength]
    plt.plot(range(len(ar)),ar)

    plt.subplot(212)
    hist(np.array(a.trace)[:,i][burnInLength::correlationLength],bins=100)    
    #plt.axvline(realFlux[i])

show()
