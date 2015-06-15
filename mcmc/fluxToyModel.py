import MCMC as mcmc
import cPickle as pickle
import numpy as np
import time
import getopt
import sys
import pandas as pd

n=4

def fluxFunction(x):
    return 1e3*x**-2.7

realFlux=[fluxFunction(i+1) for i in range(n)]

def logp(value,cov,data_observed,alpha):
    theExpectedFlux=cov.dot(np.array(value))
    log=0

    for i in range(n):
        if value[i]<0:
            return -np.inf
        var=data_observed[i]*np.log(theExpectedFlux[i]) - theExpectedFlux[i]
        log+=var

    smoothness=0

    FirstDerivative=[]
    for i in range(1,n):
        FirstDerivative.append( (np.log(value[i]) - np.log(value[i-1]))/( np.log(i+1) - np.log(i)) )

    for i in range(1,len(FirstDerivative)):
        var=alpha* np.fabs(( FirstDerivative[i] - FirstDerivative[i-1] ) / ( np.log(i+0.5) - np.log(i-0.5)))
        smoothness-=var

    return log+smoothness



# def logp(value):
#     theExpectedFlux=cov.dot(np.array(value))
#     log=0
#     for i in range(n):
#         if value[i]<0:
#             return -np.inf
#         var=data_observed[i]*np.log(theExpectedFlux[i]) - theExpectedFlux[i]
#         log+=var

#     smoothness=0

#     FirstDerivative=[]
#     for i in range(1,n):
#         FirstDerivative.append( (np.log(value[i]) - np.log(value[i-1]))/( np.log(i+1) - np.log(i)) )

#     for i in range(1,len(FirstDerivative)):
#         var=alpha* np.fabs(( FirstDerivative[i] - FirstDerivative[i-1] ) / ( np.log(i+0.5) - np.log(i-0.5)))
#         smoothness-=var

#     return log+smoothness

def defaultValues():
    cov=np.diag([0.5 for i in range(n)])
    for i in range(n):
        if i+1 < n:
            cov[i,i+1]=0.25
        if i-1 >= 0:
            cov[i,i-1]=0.25

    alpha = 0
    return cov, alpha

def generateData(cov):
    expectedAfterSmearing=cov.dot(np.array(realFlux))
    #data_observed=np.random.poisson(expectedAfterSmearing)
    data_observed=expectedAfterSmearing
    return data_observed


def main(argv):
    t0=time.time()
    filename='test.pkl'
    nStep=100000
    nThreads=4
    dirname='./'
    matrix=''

    cov, alpha = defaultValues()
    
    try:
        opts, args = getopt.getopt(argv,"f:n:t:a:d:m:")
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-f':
            filename=arg
        elif opt == '-n':
            nStep=int(arg)
        elif opt == '-t':
            nThreads=int(arg)
        elif opt == '-a':
            alpha=float(arg)
        elif opt == '-d':
            dirname=arg
        elif opt == '-m':
            # global cov
            # df=pd.read_csv('cov.txt')
            # col=list(df)
            # col.remove('Unnamed: 0')
            # offset=2
            # cov=df[col].transpose().as_matrix()[offset:n+offset,offset:n+offset]
            # print cov
            df=pd.read_csv('cov.txt')
            col=list(df)
            col.remove('Unnamed: 0')
            offset=2
            cov=df[col].transpose().as_matrix()[offset:n+offset,offset:n+offset]

    data_observed = generateData(cov)
    def wrapper(logp):
        def inner(value):
            return logp(value,cov,data_observed,alpha)
        return inner

    global logp
    logp=wrapper(logp)

    filename='alpha{}_'.format(alpha)+filename
    
    threads=[]

    for i in range(nThreads):
        theFileName=dirname+'/thread{}_'.format(i)+filename
        a=mcmc.MCMC(theFileName,initialCondition=realFlux,realValues=realFlux)
        a.setLogLikelihoodFunction(logp)
        a.setSteps(nStep)
        threads.append(a)

    for t in threads:
        print 'launching thread'
        t.start()
        time.sleep(1)

    for t in threads:
        t.join()
        
    print 'done'
    print 'time : {}'.format(time.time()-t0)

if __name__ == "__main__":
    main(sys.argv[1:])


