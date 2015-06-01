import sys
from pymc.Matplot import plot
from pymc import *
import numpy as np
from pymc import MCMC
from pylab import hist, show
import matplotlib.pyplot as plt
import pickle

#number of bins
n=4

# some response matrix
mean=[i for i in range(n)]
cov=np.diag([0.5 for i in range(n)])
for i in range(n):
    if i+1 < n:
        cov[i,i+1]=0.25
        if i-1 >= 0:
            cov[i,i-1]=0.25

#np.random.multivariate_normal(mean,cov)

def fluxFunction(x):
    return 1e3*x**-2.7


def generateMCMCFirst(dbname):
    #Generate fake data
    realFlux=[fluxFunction(i+1) for i in range(n)]
    expectedAfterSmearing=cov.dot(np.array(realFlux))
    data_observed=np.random.poisson(expectedAfterSmearing)
    print data_observed

    # Define model variables
    expected=DiscreteUniform('expectedFlux', lower=0, upper=1e4,size=n)

    

    @deterministic(plot=False)
    def expectedFluxAfterSmearing(e=expected):
        return cov.dot(np.array(e))

    observed=Poisson('observedFlux', mu=expectedFluxAfterSmearing, value=expectedAfterSmearing.astype(int), observed=True)

    model=Model({'observedFlux':observed,'expectedFlux':expected})
    M = MCMC(model, db='pickle', dbname=dbname)
    M.sample(iter=2000000, burn=100000, thin=100)
    M.db.close()
    f=open(dbname,'a')
    pickle.dump(data_observed,f)
    pickle.dump(expectedAfterSmearing,f)
    pickle.dump(realFlux,f)    


def generateMCMC(dbname):
    #Generate fake data
    realFlux=[fluxFunction(i+1) for i in range(n)]
    expectedAfterSmearing=cov.dot(np.array(realFlux))
    data_observed=np.random.poisson(expectedAfterSmearing)

    # Define model variables

    @stochastic
    def expectedFlux(value=[1000 for i in range(n)]):
        def logp(value):
            print np.array(value)
            theExpectedFlux=cov.dot(np.array(value))
            log=0
            for i in range(n):
                if value[i]<0:
                    return -np.inf
                log+=data_observed[i]*np.log10(theExpectedFlux[i]) - theExpectedFlux[i]
                

            FirstDerivative=[]
            for i in range(1,n):
                FirstDerivative.append( (np.log10(theExpectedFlux[i]) - np.log10(theExpectedFlux[i-1]))/( np.log10(i+1) - np.log10(i)) )

            smoothness=0
            alpha=0
            for i in range(1,len(FirstDerivative)):
                var=alpha* np.fabs(( np.log10(np.fabs(FirstDerivative[i])) - np.log10(np.fabs(FirstDerivative[i-1])) ) / ( np.log10(i+1) - np.log10(i)))
                smoothness-=var
            #print '{},{}'.format(log,smoothness)
                
            return log+smoothness

    model=Model({'expectedFlux':expectedFlux})
    M = MCMC(model, db='pickle', dbname=dbname)
    M.sample(iter=200000, burn=10000, thin=100)
    #M.sample(iter=20)
    M.db.close()

    f=open(dbname,'a')
    pickle.dump(data_observed,f)
    pickle.dump(expectedAfterSmearing,f)
    pickle.dump(realFlux,f)
    f.close()
    
    plot(M)
    plt.show()

generateMCMC('mcmc.pickle.tmp')

#pickle.dump(data_observed,open(dbname,'a'))

