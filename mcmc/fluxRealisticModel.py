import MCMC
from pd_model import *
import pandas as pd

def get_matrix(filename, bins):
    rres = pd.DataFrame.from_csv(filename)
    rres = rres.divide( rres.sum(axis=1), axis=0)
    rres.columns = bins
    rres.index = list(rres.index[:-1]) + [bins[-1]]
    rres.index = rres.index.astype(float)
    rres = rres.reindex(bins).fillna(0.0)
    return rres.values[:-1,:-1]


if __name__ == "__main__":
   
    rgdtBins = np.logspace(-5.0 / 19, 1, 13)
    betaBins = rgdtBins/np.sqrt(rgdtBins**2 + mp**2)

    model = PDModel(betaBins, rgdtBins)

    rgdtF = get_matrix("R_resolution.csv", rgdtBins)
    model.set_rigidity_resolution(rgdtF)
    
    betaF = get_matrix("B_resolution.csv", betaBins)
    model.set_beta_resolution(betaF)
    
    observed = np.loadtxt("observed_mock.txt")

    def logp(value):
        value = np.array(value)
        if (value < 0).any(): return -np.inf
        
        # value must be and array twice the size of binning 
        expected = model(*value.reshape((2,len(betaBins)-1)))
        log = (observed * np.log(expected) - expected).sum()

        # Didn't figure that out yet
        #firstDerivative = np.diff(np.log(value))
        #secondDerivative = np.fabs(np.diff(firstDerivative))
        #smoothness = -(alpha * secondDerivative).sum()
        
        return log #+ smoothness
   
    
