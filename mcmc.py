import sys
from pymc.Matplot import plot
from pymc import *
import numpy as np
from pymc import MCMC
from pylab import hist, show
import matplotlib.pyplot as plt

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


#Generate fake data
realFlux=[fluxFunction(i+1) for i in range(n)]
expectedAfterSmearing=cov.dot(np.array(realFlux))
data_observed=np.random.poisson(expectedAfterSmearing)

# Define model variables
expected=DiscreteUniform('expected', lower=0, upper=1e4,size=n)


@deterministic(plot=False)
def expectedFluxAfterSmearing(e=expected):
    return cov.dot(np.array(e))

observed=Poisson('observed', mu=expectedFluxAfterSmearing, value=data_observed, observed=True)

model=Model({'observedFlux':observed,'expectedFlux':expected})
M = MCMC(model)

M.sample(iter=100000, burn=10000, thin=10)

fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,10))
axes[0,0].hist(M.trace('expected')[:,0])
axes[0,1].hist(M.trace('expected')[:,1])
axes[1,0].hist(M.trace('expected')[:,2])
axes[1,1].hist(M.trace('expected')[:,3])

show()

    
