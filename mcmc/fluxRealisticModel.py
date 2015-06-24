import time
import MCMC
from pd_model import *
import pandas as pd
import numpy as np
import sys
import getopt

def load_model(filenameR, filenameB):
    frame = pd.DataFrame.from_csv(filenameR)
    frame = frame.divide(frame.sum(axis=1),axis=0)
    rgdtMeasured  = np.array(frame.columns.astype(float))
    rgdtTheoretic = np.array(frame.index.astype(float))
    rgdtF = frame.values[:-1,:-1]
    
    frame = pd.DataFrame.from_csv(filenameB)
    frame = frame.divide(frame.sum(axis=1),axis=0)
    betaMeasured  = np.array(frame.columns.astype(float))
    betaTheoretic = np.array(frame.index.astype(float))
    betaF = frame.values[:-1,:-1]

    model = PDModel(
        betaBinsTheoretic=betaTheoretic,
        betaBinsMeasured=betaMeasured ,
        rgdtBinsTheoretic=rgdtTheoretic, 
        rgdtBinsMeasured=rgdtMeasured 
        )
    
    model.set_rigidity_resolution(rgdtF)
    model.set_beta_resolution(betaF)
    return model

def make_mock_observation(fluxP, fluxD, write_file=None):
    model = load_model("datasets/R_resolution.csv", "datasets/B_resolution.csv")
    r = model(fluxP,fluxD)
    
    if not write_file is None:
        np.savetxt(write_file,r)
    
    return r

def main(argv):
    t0=time.time()
    filename='test.pkl'
    nStep=10000
    nThreads=1
    dirname='./'
    matrix=''
    alpha=0
    
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

   
    model = load_model("datasets/R_resolution.csv", "datasets/B_resolution.csv")

    #observed = np.loadtxt("datasets/observed_mock_equal.txt")

    fluxD = 100 + np.zeros(model.nBinsBT)
    fluxP = 100 + np.zeros(model.nBinsBT)
    flux = np.concatenate([fluxP,fluxD])

    observed = make_mock_observation(fluxP, fluxD)

    def logp(value):
        value = np.array(value)
        if (value < 0).any(): return -np.inf

        # value must be and array twice the size of binning 
        expected = model(*value.reshape((2,model.nBinsBT)))
        log = (observed * np.log(expected) - expected).sum()
        # Didn't figure that out yet
        #firstDerivative = np.diff(np.log(value))
        #secondDerivative = np.fabs(np.diff(firstDerivative))
        #smoothness = -(alpha * secondDerivative).sum()
        return log #+ smoothness

    sigma=1
    def proposal_function(previous_point):
        point=np.zeros(len(previous_point))
            #return previous_point+self.sigma*np.random.standard_normal(self.nVar)
        for i in range(len(previous_point)):
            while True:
                val=previous_point[i]+0.01*np.random.standard_normal()
                if val > 0:
                    point[i]=val
                    break
        return point

    filename='alpha{}_'.format(alpha)+filename
    
    threads=[]

    for i in range(nThreads):
        theFileName=dirname+'/thread{}_'.format(i)+filename
        a=MCMC.MCMC(theFileName,initialCondition=flux[:],realValues=flux[:])
        a.setProposalFunction(proposal_function)
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


