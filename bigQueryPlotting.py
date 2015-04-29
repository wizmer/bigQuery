#!/usr/bin/python
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json

bigQueryTable="[Preselected.Preselected]"

def executeQuery(theCommand):
    #os.system(theCommand)
    print theCommand
    output = subprocess.check_output(theCommand, shell=True)

    f = open('log','w')
    f.write(output)
    f.close()

    if output[0] is not '{' and output[0] is not '[':
        output=output.split('\n')[1]

    # if output[0] is '{':
    #     output = output.strip('\n')
    #     output = "[%s]" % output
    #     print output

    try:
        jsonData=json.loads(output)

    except ValueError:
        print 'no json data'


def hist( nBins, firstBin, lastBin, var, cut='' ):
    binWidth=float(lastBin-firstBin)/nBins
    theCommand="bq --format json query -n " + str(nBins+10) + " \'SELECT INTEGER(FLOOR( ( (" + var + ") - (" + str(firstBin) + "))/(" + str(binWidth) + ") )) as binX,COUNT(1) FROM " + bigQueryTable
    if cut :
        theCommand+=" WHERE (" + cut + ") "

    theCommand+=" GROUP BY binX"
    theCommand+=" HAVING binX >= 0 AND binX < (" + str(nBins)+ " )"
    theCommand+=" ORDER BY binX \'"

    try:
        data = executeQuery(theCommand)

        L=[0 for i in range(nBins)]


        for d in data:
            L[int(d['binX'])]=float(d['f0_'])
            
            x=[firstBin + i*binWidth for i in range(nBins)]
            hist, bins = np.histogram(x, range=(firstBin,lastBin),bins=nBins, weights=L)
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            #    plt.yscale('log')
            plt.bar(center, hist, width=width)
            plt.show()
        except Exception:
            print 'no json data'
                


def hist2d( nBinsX, firstBinX, lastBinX, nBinsY, firstBinY, lastBinY, varX, varY, cut='' ):
    binWidthX=float(lastBinX-firstBinX)/nBinsX
    binWidthY=float(lastBinY-firstBinY)/nBinsY

    theCommand="bq --format json query -n " + str(nBinsX*nBinsY+10) + " \'SELECT INTEGER(FLOOR( ((" + varX + ") - (" + str(firstBinX) + "))/(" + str(binWidthX) + ") )) as binX, INTEGER(FLOOR( ((" + varY + ") - (" + str(firstBinY) + "))/(" + str(binWidthY) + "))) as binY,COUNT(1) FROM " + bigQueryTable
    if cut :
        theCommand+=" WHERE (" + cut + ")"

    theCommand+=" GROUP BY binX, binY"
    theCommand+=" HAVING binX >= 0 AND binX < (" + str(nBinsX) + ") AND binY >= 0 AND binY < (" + str(nBinsY) + ")"
    theCommand+=" ORDER BY binX,binY \'"


    L=[0 for i in range(nBinsX*nBinsY)]

    data = executeQuery(theCommand)
    
    for d in data:
        L[int(d['binX'])*nBinsY+int(d['binY'])]=float(d['f0_'])
        
    x=[firstBinX + (i/nBinsY)*binWidthX  for i in range(nBinsX*nBinsY)]
    y=[firstBinY + (i%nBinsY)*binWidthY for i in range(nBinsX*nBinsY)]

    plt.hist2d(x,y, bins=[nBinsX,nBinsY], range=[[firstBinX, lastBinX], [firstBinY, lastBinY]], weights=L)
    plt.colorbar()

    plt.show()


def setTable(newTableName):
    global bigQueryTable
    bigQueryTable=newTableName


def makeSelectionMask(cutList):
    data=executeQuery("bq --format json show Preselected.Preselected")
    for d in data['schema']['fields']:
        if d['name'] == 'selStatus':
            selStatusName=d['description']

    selStatusName=selStatusName.split(',')
    bitIndex=dict()
    for i in range(len(selStatusName)):
        bitIndex[selStatusName[i]]=i

    theMask=0
    for cut in cutList:
        if cut in bitIndex:
            theMask+=1<<bitIndex[cut]
        else:
            print 'Cut not found !'
            return -1

    return theMask
