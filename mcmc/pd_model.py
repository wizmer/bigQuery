import numpy as np

__all__ = ["PDModel", "mp", "md"]

# Masses
mp = 0.9382
md = 1.8756 


def R_from_beta(beta, m): return m*beta/np.sqrt(1-beta*beta)
def beta_from_R(R, m): return R/np.sqrt(R*R+m*m) 


import itertools
def iterpairs(l):
    i1,i2 = itertools.tee(l)
    next(i2)
    return((x1,x2) for x1,x2 in itertools.izip(i1,i2))

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def _delta_matrix(delta, betaBins, rgdtBins, m):
    cells = itertools.product(
        enumerate(iterpairs(betaBins)), 
        enumerate(iterpairs(rgdtBins))
    )
    for (nb, (bmin, bmax)), (nR, (Rmin, Rmax)) in cells:
        b1,b2 = beta_from_R(Rmin, m), beta_from_R(Rmax, m)
        ovlp = getOverlap((bmin,bmax), (b1,b2))/(bmax-bmin)
        delta[nb,nR] = ovlp

def _test_matrix_shape(matrix, (size1, size2), name):
    if matrix.ndim != 2:
        raise ValueError(name + ' must be 2D')
    if matrix.shape[0] != size1:
        raise ValueError(name + ' first dimension mismatch')
    if matrix.shape[1] != size2: 
        raise ValueError(name + ' second dimension mismatch')
        
        
class PDModel(object):
    def __init__(self, betaBinsTheoretic = None, 
                       betaBinsMeasured  = None, 
                       rgdtBinsTheoretic = None,
                       rgdtBinsMeasured  = None):

        if betaBinsTheoretic is None: raise ValueError('Missing betaBinsTheoretic ')
        if betaBinsMeasured  is None: raise ValueError('Missing betaBinsMeasured  ')
        if rgdtBinsTheoretic is None: raise ValueError('Missing rgdtBinsTheoretic ')
        if rgdtBinsMeasured  is None: raise ValueError('Missing rgdtBinsMeasured  ')
        
        if betaBinsTheoretic.ndim != 1: raise ValueError('betaBinsTheoretic must be 1D')
        if betaBinsMeasured .ndim != 1: raise ValueError('betaBinsMeasured  must be 1D')
        if rgdtBinsTheoretic.ndim != 1: raise ValueError('rgdtBinsTheoretic must be 1D')
        if rgdtBinsMeasured .ndim != 1: raise ValueError('rgdtBinsMeasured  must be 1D')
        
        self.nBinsBT   = betaBinsTheoretic.shape[0] - 1
        self.betaBinsT = betaBinsTheoretic.copy()
        self.nBinsBM   = betaBinsMeasured.shape[0] - 1
        self.betaBinsM = betaBinsMeasured.copy()
        self.betaF     = None
        
        self.nBinsRT   = rgdtBinsTheoretic.shape[0] - 1 
        self.rgdtBinsT = rgdtBinsTheoretic.copy() 
        self.nBinsRM   = rgdtBinsMeasured.shape[0] - 1 
        self.rgdtBinsM = rgdtBinsMeasured.copy() 
        self.rgdtF     = None 
        
        self.deltaP = np.zeros((self.nBinsBT, self.nBinsRT))
        self.deltaD = np.zeros((self.nBinsBT, self.nBinsRT))
        
        _delta_matrix(self.deltaPT, self.betaBinsT, self.rgdtBinsT, mp)
        _delta_matrix(self.deltaDT, self.betaBinsT, self.rgdtBinsT, md)
        
    
    def set_rigidity_resolution(self, rgdtF):
        _test_matrix_shape(rgdtF, [self.nBinsRT, self.nBinsRM], 'Rigidity resolution matrix')
        self.rgdtF = rgdtF.copy()
        
    def set_beta_resolution(self, betaF):
        _test_matrix_shape(betaF, [self.nBinsBT, self.nBinsBM], 'Beta resolution matrix')
        self.betaF = betaF.copy()

#    Exposure time not ready yet   
#    def set_exposure_time(self, latitudes, exptime):
#        if latitudes.ndim != 1:
#            raise ValueError('Latitude bins dimension must be 1')
#        self.nBinsLat = latitudes.shape[0] - 1
#        
#        _test_matrix_shape(exptime, (self.nBinsLat, self.nBinsR), 'Exposure time matrix')
#        self.exptime = exptime.copy()
        
    def __call__(self, fluxP, fluxD):
        if fluxP.ndim != 1:
            raise ValueError('Proton flux must be 1D')
        if fluxP.shape[0] != self.nBinsBT:
            raise ValueError('Proton flux size doesnt correspond to beta bins')
        if fluxD.ndim != 1:
            raise ValueError('Deuteron flux must be 1D')
        if fluxD.shape[0] != self.nBinsBT:
            raise ValueError('Deuteron flux size doesnt correspond to beta bins')
       
        # Making 2D(beta,R)
        fluxP = self.deltaP * fluxP[:,np.newaxis]
        fluxD = self.deltaD * fluxD[:,np.newaxis]
        
#       # Making 3D(lat, beta, R)
#       fluxP = fluxP[np.newaxis,:,:] * self.exptime[:,np.newaxis,:] 
#       fluxD = fluxD[np.newaxis,:,:] * self.exptime[:,np.newaxis,:] 
       
        # Smearing beta (note that tensordot transposes dimensions)
        fluxP = np.tensordot(fluxP, self.betaF, axes=([0],[0]))
        fluxD = np.tensordot(fluxD, self.betaF, axes=([0],[0]))
        # Smearing Rself
        fluxP = np.tensordot(fluxP, self.rgdtF, axes=([0],[0]))
        fluxD = np.tensordot(fluxD, self.rgdtF, axes=([0],[0]))
        
        return fluxP+fluxD
        
        

        
        
        
        
        
        
        
