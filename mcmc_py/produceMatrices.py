import pandas as pd
from pd_model import *
import seaborn as s
from matplotlib.colors import LogNorm
import bq
from histQueryFactory import *
import matplotlib.pyplot as plt
import numpy as np

selectionStatus ="1972401"
#selectionStatus="2097137"
#selectionStatus=2097151

s.set(rc={'image.cmap': "jet"})
#matplotlib.figsize(12,10)
# rcParams['figure.facecolor'] = (1,1,1,1)
# rcParams['savefig.facecolor'] = (1,1,1,1)

################################################################################
##
##  MAKE BINNINGS
##
################################################################################
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

#rcParams['savefig.dpi'] = 160
#plt.figsize(8,4)
# for b in betaTheoretic:plot([0,10],[b,b],'k',lw=0.5)
# for R in rgdtTheoretic:plot([R,R], [0.5,1],'k',lw=0.5)
# xlabel("R theoretic")
# ylabel("$\\beta$ theoretic")
# xlim(rgdtTheoretic[0],10)
# ylim(0.49,1)

# x = np.linspace(0,10,100)
# plot(x, beta_from_R(x,mp),'r')
# plot(x, beta_from_R(x,md),'b')

#
# Measured values binnings are not very fancy
#
betaMeasured = 1/np.linspace(0.5,2,28)
betaMeasured.sort()

rgdtMeasured = np.logspace(-5.0 / 19, 1, 25)

# rcParams['savefig.dpi'] = 160
#plt.figsize(8,4)
# for b in betaMeasured:plot([0,10],[b,b],'k',lw=0.5)
# for R in rgdtMeasured:plot([R,R], [0.5,1.2],'k',lw=0.5)
# xlabel("R measured")
# ylabel("$\\beta$ measured")
# xlim(rgdtTheoretic[0],10)
# ylim(0.49,1.21)

# x = np.linspace(0,10,100)
# plot(x, beta_from_R(x,mp),'r')
# plot(x, beta_from_R(x,md),'b')

################################################################################
##
##  MAKING BQ CLEAR
##
################################################################################

client = bq.Client.Get()
schema = client.GetTableSchema({
        'projectId': 'ams-test-kostya',
        'datasetId': 'AMS',
        'tableId': 'protonsB1034'
    })

bitFields = None
for field in schema['fields']:
    if field['name'] != 'selStatus':
        continue
    bitFields = field['description'].split(',')
    break

get_cumulative_mask = lambda sel: (1 << (bitFields.index(sel))) - 1
                            

import itertools
def iterpairs(l):
    i1,i2 = itertools.tee(l); next(i2)
    return((x1,x2) for x1,x2 in itertools.izip(i1,i2))

def build_case_string(var, name, bins):
    casestr = ["CASE"]
    casestr += ["  WHEN {1} < {0} AND {0} < {2} THEN {3}".format(var,s,e,n)
                for n,(s,e)  in enumerate(iterpairs(bins))]
    casestr += ["  ELSE -1",  "END as {0}".format(name)]
    return '\n'.join(casestr)

################################################################################
##
##  BETA MATRIX
##
################################################################################
vs =  ",\n".join([ build_case_string("BetaTOF", "B_bin", betaMeasured),
                   build_case_string("GenMomentum/SQRT(0.88022 + POW(GenMomentum,2))", "Gen_bin", betaTheoretic),
                   "COUNT(1) as count" ])

h = "SELECT\n" + vs + """
FROM
   AMS.protonsB1034
WHERE
   selStatus&""" + selectionStatus + """=""" + selectionStatus + """
GROUP BY B_bin,Gen_bin
ORDER BY B_bin,Gen_bin"""

tableid = client.Query(str(h))['configuration']['query']['destinationTable']
bq_table = client.ReadTableRows(tableid)

frame = pd.DataFrame(bq_table, columns=['Bbin', 'GenBin', 'Count']).astype(int)
frame['Bbin'] = frame['Bbin'].map(lambda x: betaMeasured[x] )
frame['GenBin'] = frame['GenBin'].map(lambda x: betaTheoretic[x] )
frame = frame.set_index(list(frame.columns[:-1])).unstack()['Count'].fillna(0)

frame.T.to_csv("./datasets/B_resolution.csv")

#plt.figsize(5,4)
# cc = plot_matrix(frame.ix[frame.index[:-1]][frame.columns[:-1]],
#                             norm=LogNorm(vmin=1,vmax=10**6))
# colorbar(cc)
# xlabel("$\\beta$ True")
# ylabel("$\\beta$ Measured")
# ylim(0.6,1.4)

################################################################################
##
##  RIGIDITY MATRIX
##
################################################################################
vs =  ",\n".join([ build_case_string("R", "R_bin", rgdtMeasured),
                                      build_case_string("GenMomentum", "Gen_bin", rgdtTheoretic),
                                      "COUNT(1) as count" ])

h = "SELECT\n" + vs + """
FROM
   AMS.protonsB1034
WHERE
   selStatus&""" + selectionStatus + """=""" + selectionStatus + """
GROUP BY R_bin,Gen_bin
ORDER BY R_bin,Gen_bin"""

tableid = client.Query(str(h))['configuration']['query']['destinationTable']
bq_table = client.ReadTableRows(tableid)
frame = pd.DataFrame(bq_table, columns=['Rbin', 'GenBin', 'Count']).astype(int)
frame['Rbin'] = frame['Rbin'].map(lambda x: rgdtMeasured[x] )
frame['GenBin'] = frame['GenBin'].map(lambda x: rgdtTheoretic[x] )
frame = frame.set_index(list(frame.columns[:-1])).unstack()['Count'].fillna(0)
frame.T.to_csv("./datasets/R_resolution.csv")
#plt.figsize(5,4)
# cc = plot_matrix(frame.ix[frame.index[:-1]][frame.columns[:-1]],
#                             norm=LogNorm(vmin=1,vmax=10**6))
# colorbar(cc)
# xlabel("$R$ True")
# ylabel("$R$ Measured")
# xlim(0.5,10)

################################################################################
##
##  TARGET
##
################################################################################

vs =  ",\n".join([ build_case_string("R", "R_bin", rgdtMeasured),
                                      build_case_string("BetaTOF", "B_bin", betaMeasured),
                                      "COUNT(1) as count" ])

h = "SELECT\n" + vs + """
FROM
   AMS.Data
WHERE
   selStatus&""" + selectionStatus + """=""" + selectionStatus + """ AND
   Livetime > 0.5 AND
   ABS(Latitude) > 1.0
GROUP BY R_bin,B_bin
ORDER BY R_bin,B_bin"""
tableid = client.Query(str(h))['configuration']['query']['destinationTable']
bq_table = client.ReadTableRows(tableid)

frame = pd.DataFrame(bq_table, columns=['Rbin', 'Beta', 'Count']).astype(int)
frame['Rbin'] = frame['Rbin'].map(lambda x: rgdtMeasured[x] )
frame['Beta'] = frame['Beta'].map(lambda x: betaMeasured[x] )
frame = frame.set_index(list(frame.columns[:-1])).unstack()['Count'].fillna(0)
frame = frame.T

#plt.figsize(13,6)
# cc = plot_matrix(frame.ix[frame.index[:-1]][frame.columns[:-1]],
#                             norm=LogNorm(vmin=1,vmax=10**6))
# colorbar(cc)
# xlabel("$R$ Measured")
# ylabel("$\\beta$ Measured")
# xlim(0.5,5)
# ylim(0.5,1.3)

np.savetxt("datasets/observed_data_selStatus"+selectionStatus+".txt",frame.values)
