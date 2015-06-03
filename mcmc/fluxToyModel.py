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
        a.setSteps(100000)
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
    f=open('thread0_'+filename,'r')
    metadata=pickle.load(f)
    f.close()
    # realFlux=pickle.load(f)
    # data_observed=pickle.load(f)
    print metadata

    burnInLength=0
    correlationLength=1

    totalCount=np.zeros((4,100))
    heatmap=np.zeros((4,4,100,100))
    extent=np.zeros((4,4,4))

    bins=np.zeros((4,101))
    firstHisto=True

    for iFile in range(4):
        f=open('thread{}_'.format(iFile)+filename,'r')
        pickle.load(f) #skip metadata
        remainingBurnInLength=burnInLength
        while True:
            try:
                chunk=pickle.load(f)
                print 'chunk loaded'

                # skip chunk before burn-in length
                if remainingBurnInLength >= chunk['numberOfSteps']:
                    print 'skip chunk'
                    remainingBurnInLength-=chunk['numberOfSteps']
                    if remainingBurnInLength < 0: remainingBurnInLength=0
                    continue
            except:
                break

            try:
                if firstHisto:
                    for iVar in range(metadata['nVar']):
                        counts, bins[iVar], totalPatch=hist(np.array(chunk['trace'])[:,iVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],bins=100)
                        totalCount[iVar]+=counts
                        if iVar == 0:
                            print counts

                        for jVar in range(metadata['nVar']):
                            heatmap[iVar][jVar], xedges, yedges = np.histogram2d(np.array(chunk['trace'])[:,iVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],np.array(chunk['trace'])[:,jVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],bins=100)
                            extent[iVar][jVar] = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

                    firstHisto=False

                else:
                    for iVar in range(metadata['nVar']):
                        counts, b, totalPatch=hist(np.array(chunk['trace'])[:,iVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],bins=bins[iVar])
                        totalCount[iVar]+=counts
                        if iVar == 0:
                            print counts

                        for jVar in range(metadata['nVar']):
                            h, xedges, yedges = np.histogram2d(np.array(chunk['trace'])[:,iVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],np.array(chunk['trace'])[:,jVar][remainingBurnInLength:chunk['numberOfSteps']:correlationLength],bins=100,range=[[extent[iVar][jVar][0],extent[iVar][jVar][1]],[extent[iVar][jVar][2],extent[iVar][jVar][3]]])
                            heatmap[iVar][jVar]+=h

            except:
                pass

            remainingBurnInLength-=chunk['numberOfSteps']
            if remainingBurnInLength < 0: remainingBurnInLength=0

            del chunk

    plt.figure('Correlation plots')
    for iVar in range(metadata['nVar']):
        for jVar in range(metadata['nVar']):
            plotNum=iVar*metadata['nVar']+jVar+1
            plt.subplot(metadata['nVar'],metadata['nVar'],plotNum)
            if iVar==jVar:
                plt.plot(bins[iVar][:-1],totalCount[iVar],drawstyle='steps')
                plt.axvline(realFlux[iVar])
            elif iVar>jVar:
                plt.imshow(heatmap[jVar][iVar], extent=extent[jVar][iVar],aspect="auto")
                plt.title('titre {},{}'.format(iVar,jVar))
    
    show()

    return mcmc

#a=generateMCMC('test.pkl',4)
a=loadMCMC('test.pkl')


