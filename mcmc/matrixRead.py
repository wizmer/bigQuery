#
# "Theoretical" binning (the one we will look the spectrum for)
#  is bootstrapped from R-vs-beta curves for D and P 
#
def make_beta_bins(beta):
    bbins = []
    for i in range(10):
        bbins.append(beta)
        beta = beta_from_R(R_from_beta(beta,md),mp)
    return bbins

bbins = make_beta_bins(0.5)
bbins += make_beta_bins((bbins[1]+bbins[0])/2)
bbins = sorted(bbins)

mid1,mid2 = (bbins[1]+bbins[0])/2,(bbins[2]+bbins[1])/2
bbins += make_beta_bins(mid1)
bbins += make_beta_bins(mid2)
bbins = np.array(sorted(bbins))

betaTheoretic, rgdtTheoretic = np.array([bbins, R_from_beta(bbins, mp)])



#
# Measured values binnings are not very fancy
#
betaMeasured = 1/np.linspace(0.5,2,28)
betaMeasured.sort()

rgdtMeasured = np.logspace(-5.0 / 19, 1, 25)

# This reads matrices

def get_matrix(filename, binsT, binsM):
    rres = pd.DataFrame.from_csv(filename)
    rres = rres.divide( rres.sum(axis=1), axis=0)
    rres.columns = binsM
    rres.index = list(rres.index[:-1]) + [binsT[-1]]
    rres.index = rres.index.astype(float)
    rres = rres.reindex(binsT).fillna(0.0)
    return rres.values[:-1,:-1]

rgdtF = get_matrix("datasets/R_resolution.csv", rgdtTheoretic, rgdtMeasured)
model.set_rigidity_resolution(rgdtF)

betaF = get_matrix("datasets/B_resolution.csv", betaTheoretic, betaMeasured)
model.set_beta_resolution(betaF)
