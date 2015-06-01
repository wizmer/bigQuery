import sys
from pymc import *
import numpy as np
from pymc import MCMC
from pylab import hist, show
import matplotlib.pyplot as plt
import pickle

#number of bins
n=10

# def plot(M):
#     # fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,10))
#     # axes[0,0].hist(M.trace('expected')[:,0])
#     # axes[0,1].hist(M.trace('expected')[:,1])
#     # axes[1,0].hist(M.trace('expected')[:,2])
#     # axes[1,1].hist(M.trace('expected')[:,3])
#     for i in range(n):
#         plt.figure()
#         hist(M.trace('expected')[:,i])


def loadDB(filename):
    unpickled = []

    db = pymc.database.pickle.load(filename)
    f=open(filename)
    pickle.load(f)
    data_observed=pickle.load(f)
    expectedFluxAfterSmearing=pickle.load(f)
    realFlux=pickle.load(f)

    # print realFlux
    # print data_observed
    # print expectedFluxAfterSmearing

    for i in range(n):
        plt.figure()
        hist(db.trace('expectedFlux')[:,i],bins=100)
        plt.axvline(realFlux[i])

    print realFlux

    plt.figure()
    plt.scatter(x=db.trace('expectedFlux')[:,0],y=db.trace('expectedFlux')[:,1])

    show()

loadDB('mcmc.pickle.tmp')


