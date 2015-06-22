import MCMC
from pd_model import *
import pandas as pd

if __name__ == "__main__":
    
    frame = pd.DataFrame.from_csv("datasets/R_resolution.csv")
    frame = frame.divide(frame.sum(axis=1),axis=0)
    rgdtMeasured  = np.array(frame.columns.astype(float))
    rgdtTheoretic = np.array(frame.index.astype(float))
    rgdtF = frame.values[:-1,:-1]

    frame = pd.DataFrame.from_csv("datasets/B_resolution.csv")
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

    observed = np.loadtxt("datasets/observed_mock_equal.txt")

    def logp(value):
        value = np.array(value)
        if (value < 0).any(): return -np.inf
        
        # value must be and array twice the size of binning 
        expected = model(*value.reshape((2,len(betaTheoretc)-1)))
        log = (observed * np.log(expected) - expected).sum()

        # Didn't figure that out yet
        #firstDerivative = np.diff(np.log(value))
        #secondDerivative = np.fabs(np.diff(firstDerivative))
        #smoothness = -(alpha * secondDerivative).sum()
        
        return log #+ smoothness
   
    flux = 1 + np.zeros_like(betaTheoretic)[:-1]
    flux = np.concatenate([flux,flux])
    
    mcmc = MCMC.MCMC("testtest", flux[:], flux[:]) 
    mcmc.setLogLikelihoodFunction(logp)
    mcmc.setSteps(10000) 

    mcmc.start()
    mcmc.join()