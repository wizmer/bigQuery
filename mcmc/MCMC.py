import pickle
import numpy as np
import matplotlib.pyplot as plt
import os

n=4


mean=[1,400,30,50]

def default_logp(value):
    log=0
    for i in range(n):
        var=-((value[i]-mean[i])/mean[i])**2
        log+=var
        return log
    
class MCMC:
    def __init__(self,filename):
        self.filename=filename
        self.trace=[]
        self.log_likelihood=[] # the chain containing the log likelihood of all accepted points
        self.nStep=10
        self.current_point=[10 for i in range(n)] # the chain starting point
        self.nVar=n
        self.verbose=False
        self.sigma=[1 for i in range(n)] # the sigmas for the proposal function
        self.logp=default_logp
        self.chunk=[self.log_likelihood,self.trace]

    def saveMetaData(self):
        print 'saving metadata...'
        f=open(self.filename+'.pkl','a')
        pickle.dump(self.metadata,f)
        f.close()

    def saveChunk(self):
        print 'saving chunk...'
        f=open(self.filename+'.pkl','a')
        pickle.dump(self.chunk,f)
        f.close()
        del self.log_likelihood[:]
        del self.trace[:]

    def setInitialConditions(self,initialCondition):
        self.current_point=initialCondition

    def setSigmas(self, sigma):
        self.sigma=sigma

    def proposal_function(self,previous_point):
        cov=np.diag(self.sigma)
        return np.random.multivariate_normal(previous_point,cov)

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
        self.trace.append(self.current_point)
        self.log_likelihood.append(self.current_log_likelihood)

    def loop(self):
        self.metadata={'initialPoints':self.current_point,'nStep':self.nStep,'sigmas':self.sigma,'nVar':self.nVar}
        os.system('rm -f '+filename+'.pkl')
        self.saveMetaData()

        self.current_log_likelihood=self.logp(self.current_point)

        for i in range(self.nStep):
            if i%50000 == 0: print i
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

            if len(self.log_likelihood) >= 5000000:
                self.saveChunk()

        if len(self.log_likelihood) > 0: self.saveChunk()


def generateMCMC(filename):
    a=MCMC()
    a.setSteps(1000000)
    a.loop()
    f=open(filename,'wb')
    pickle.dump(a,f)
    f.close()
    return a

def loadMCMC(filename):
    f=open(filename,'rb')
    return pickle.load(f)
