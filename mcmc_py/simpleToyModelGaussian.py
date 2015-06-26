import MCMC as mcmc
import pickle
import matplotlib.pyplot as plt
import numpy as np
from pylab import hist, show

n=4
mean=[1,400,30,50]
sigma=[3,20,5,5]


def logp(value):
    log=0
    for i in range(n):
        var=-((value[i]-mean[i])/sigma[i])**2
        log+=var
    return log

def generateMCMC(filename):
    a=mcmc.MCMC()
    a.setSteps(100000)
    a.setSigmas(sigma)
    a.current_point=mean
    a.setLogLikelihoodFunction(logp)
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

burnInLength=0
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
