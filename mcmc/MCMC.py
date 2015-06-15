import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import os
import gzip
from threading import Thread
from multiprocessing import Process
import pandas as pd
import time

mean=[1,400,30,50]

def default_logp(value):
    log=0
    for i in range(4):
        var=-((value[i]-mean[i])/mean[i])**2
        log+=var
        return log
    
class MCMC(Process):
    chunkSize=1e6
    
    def __init__(self,filename,initialCondition,realValues=[]):
        Process.__init__(self)
        self.realValues=realValues
        self.filename=filename
        self.nStep=10
        self.current_point=initialCondition # the chain starting point
        self.nVar=len(initialCondition)
        self.verbose=False
        self.trace=np.zeros((self.chunkSize,self.nVar))
        self.log_likelihood=np.zeros(self.chunkSize) # the chain containing the log likelihood of all accepted points
        self.sigma=[1 for i in range(self.nVar)] # the sigmas for the proposal function
        self.logp=default_logp
        self.chunkStepNumber=0
        self.seed=0

    def run(self):
        self.loop()

    def saveMetaData(self):
        print 'save metadata'
        f=open(self.filename,'wb')
        pickle.dump(self.metadata,f,2)
        f.close()

    def saveChunk(self):
        print 'saving chunk...'
        #f=gzip.open(self.filename,'a')
        df=pd.DataFrame(self.trace[:self.chunkStepNumber])
        dfLog=pd.DataFrame(self.log_likelihood[:self.chunkStepNumber],columns=['log'])
        dfTotal=df.join(dfLog)
        print 'df done'
        f=open(self.filename,'ab')
        print 'file opened'

        pickle.dump(dfTotal, f, protocol=2)
          
        f.close()
        print 'file closed'
        self.log_likelihood=np.zeros(self.chunkSize)
        self.trace=np.zeros((self.chunkSize,self.nVar))
        self.chunkStepNumber=0
        print 'chunk saved'

    def setSigmas(self, sigma):
        self.sigma=sigma

    def proposal_function(self,previous_point):
        return previous_point+self.sigma*np.random.standard_normal(self.nVar)

    def setProposalFunction(self,f):
        self.proposal_function=f

    def setLogLikelihoodFunction(self,f):
        self.logp=f

    def setSteps(self,steps):
        self.nStep=steps

    def setVerbose(self,isVerbose):
        self.verbose=isVerbose

    def updateStep(self):
        if self.verbose: print 'accepted'
        self.current_point=self.proposed_point
        self.current_log_likelihood=self.proposed_log_likelihood
        self.trace[self.chunkStepNumber]=self.current_point
        self.log_likelihood[self.chunkStepNumber]=self.current_log_likelihood
        self.chunkStepNumber+=1

    def loop(self):
        np.random.seed()

        self.metadata={'initialPoints':self.current_point,'nStep':self.nStep,'sigmas':self.sigma,'nVar':self.nVar,'realValues':self.realValues}
        self.saveMetaData()

        self.current_log_likelihood=self.logp(self.current_point)

        for i in range(self.nStep):
            if i%100000 == 0: print '{}/{} : {}%'.format(i,self.nStep,int(float(i)/float(self.nStep)*100))
            self.proposed_point=self.proposal_function(self.current_point)
            self.proposed_log_likelihood=self.logp(self.proposed_point)
            
            self.the_likelihood_ratio=np.exp(self.proposed_log_likelihood-self.current_log_likelihood)

            if self.verbose:
                print 'current: {}'.format(self.current_point)
                print 'proposed: {}'.format(self.proposed_point)
                
            if self.the_likelihood_ratio > 1:
                self.updateStep()
                
            else:
                if np.random.rand() < self.the_likelihood_ratio:
                    self.updateStep()

            if self.verbose:
                print 'current_log_likelihood: {}'.format(self.current_log_likelihood)    
                print 'proposed_log_likelihood: {}'.format(self.proposed_log_likelihood)
                print 'likelihood_ratio: {}'.format(self.the_likelihood_ratio)

            if self.chunkStepNumber >= self.chunkSize:
                self.saveChunk()


        print 'step number'+str(self.chunkStepNumber)
        if self.chunkStepNumber > 0: self.saveChunk()


