import MCMC as mcmc
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from pylab import hist, show
import gzip
import time

n=4

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
#data_observed=np.random.poisson(expectedAfterSmearing)
data_observed=expectedAfterSmearing

def logp(value):
    theExpectedFlux=cov.dot(np.array(value))
    log=0
    for i in range(n):
        if value[i]<0:
            return -np.inf
        var=data_observed[i]*np.log(theExpectedFlux[i]) - theExpectedFlux[i]
        log+=var

    smoothness=0
    alpha=0

    # FirstDerivative=[]
    # for i in range(1,n):
    #     FirstDerivative.append( (np.log(theExpectedFlux[i]) - np.log(theExpectedFlux[i-1]))/( np.log(i+1) - np.log(i)) )

    # for i in range(1,len(FirstDerivative)):
    #     var=alpha* np.fabs(( FirstDerivative[i] - FirstDerivative[i-1] ) / ( np.log(i+0.5) - np.log(i-0.5)))
    #     smoothness-=var

#    print '{},{}'.format(log,smoothness)

    return log+smoothness

def generateMCMC(filename,threadNumber=1):
    threads=[]
    for i in range(threadNumber):
        a=mcmc.MCMC('thread{}_'.format(i)+filename)
        a.setLogLikelihoodFunction(logp)
        a.setInitialConditions(realFlux)
        a.setSteps(1000000)
        threads.append(a)

    for t in threads:
        print 'launching thread'
        t.start()
        time.sleep(1)

    for t in threads:
        t.join()
        
    print 'done'

    # f=open(filename,'wb')
    # pickle.dump(a,f)
    # pickle.dump(realFlux,f)
    # pickle.dump(data_observed,f)
    # f.close()
    return a

def loadMCMC(filename):
    f=open(filename,'r')
    metadata=pickle.load(f)
    # realFlux=pickle.load(f)
    # data_observed=pickle.load(f)
    print metadata

    burnInLength=0
    correlationLength=1

    totalCount=np.zeros((4,100))
    bins=np.zeros((4,101))
    firstHisto=[True]*metadata['nVar']
    
    while True:
        try:
            chunk=pickle.load(f)
            print 'chunk loaded'
            if burnInLength >= chunk['numberOfSteps']:
                print 'skip chunk'
                burnInLength-=chunk['numberOfSteps']
                if burnInLength < 0: burnInLength=0

                continue
        except:
            break
        #hist(np.array(a.trace)[:,i][burnInLength::correlationLength],bins=100)

        try:
            for iVar in range(metadata['nVar']):
                if firstHisto[iVar]:
                    counts, bins[iVar], totalPatch=hist(np.array(chunk['trace'])[:,iVar][burnInLength:chunk['numberOfSteps']:correlationLength],bins=100)
                    totalCount[iVar]+=counts
                    firstHisto[iVar]=False
                else:
                    counts, b, totalPatch=hist(np.array(chunk['trace'])[:,iVar][burnInLength:chunk['numberOfSteps']:correlationLength],bins=bins[iVar])
                    totalCount[iVar]+=counts
        except:
            pass


        burnInLength-=chunk['numberOfSteps']
        if burnInLength < 0: burnInLength=0

        del chunk


    for iVar in range(metadata['nVar']):
        plt.subplot(2,2,iVar+1)
        plt.plot(bins[iVar][:-1],totalCount[iVar],drawstyle='steps')
        plt.axvline(realFlux[iVar])
    show()

    return mcmc

#a=generateMCMC('test.pkl',4)
a=loadMCMC('thread0_test.pkl')


# for i in range(n):
#     plt.figure('Bin #{}'.format(i))
#     plt.subplot(211)

#     ar=np.array(a.trace)[:,i][burnInLength::correlationLength]
#     plt.plot(range(len(ar)),ar)
#     plt.axhline(realFlux[i])

#     plt.subplot(212)
#     hist(np.array(a.trace)[:,i][burnInLength::correlationLength],bins=100)    
#     plt.axvline(realFlux[i])

# plt.figure('Correlation plots')
# for i in range(n):
#     for j in range(n):
#         plotNum=i*n+j+1
#         plt.subplot(n,n,plotNum)
#         if i==j:
#             hist(np.array(a.trace)[:,i][burnInLength::correlationLength],bins=100)
#             plt.axvline(realFlux[i])
#         elif i>j:
#             plt.plot(np.array(a.trace)[:,j][burnInLength::correlationLength],np.array(a.trace)[:,i][burnInLength::correlationLength])
#             plt.title('titre {},{}'.format(i,j))

# plt.figure('Log likelihood')
# ar=a.log_likelihood[burnInLength::correlationLength]
# plt.plot(range(len(ar)),ar)

# plt.figure('flux')
# plt.subplot(211)
# ar=np.array(a.trace)[burnInLength::correlationLength]
# for i in range(len(ar)):
#     plt.plot(range(len(ar[i])),ar[i])

# plt.subplot(212)
# plt.xscale('log')
# plt.yscale('log')
# for i in range(len(ar)):
#     plt.plot(range(1,len(ar[i])+1),ar[i])
# show()
