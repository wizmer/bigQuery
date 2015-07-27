#!/usr/bin/python
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json
import bq
import pickle as pkl
import histQueryFactory
import pandas as pd
from matplotlib import cm
from matplotlib.pylab import *
from scipy.optimize import curve_fit

bigQueryTable="AMS.Data"

def executeQuery(theCommand):
    #os.system(theCommand)
    print theCommand
    try:

        rawOutput = subprocess.check_output(theCommand, stderr=subprocess.STDOUT,shell=True)

        f = open('log','w')
        f.write(rawOutput)
        f.close()
        
        output=rawOutput
        if output[0] is not '{' and output[0] is not '[':
            output=rawOutput.split('\n')[1]

        jsonData=json.loads(output)
        print rawOutput[:1000]
        return jsonData

    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode)
        outputError=exc.output.split('\n')
        for o in outputError:
            print o

        if exc.returncode == 2 and theCommand.find('--require_cache') != -1 :
            print 'Table not found'
            print 'Do you want to retry without the \'--require_cache\' option ?'
            inp=raw_input('[Y,n] ? ')
            if inp!='Y' and inp!='y' and inp!='':
                print 'Doing nothing then'
            else:
                return executeQuery(theCommand.replace('--require_cache',''))

        return dict()
        
    except ValueError:
        print 'no json data'
        print rawOutput
        return json.loads('[{"f0_":"0","binX":"0","binY":"0"}]')
    return dict()

      
# def executeQueryNew(theCommand):
#     client = bq.Client.Get()
#     tableid = client.Query(theCommand)['configuration']['query']['destinationTable']
#     table = client.ReadTableRows(tableid)
#     print table
#     return table

def binCenterFromArray(var,array,aliasName='binX'):
    res='CASE\n'
    for i in range(len(array)-1):
        res+='    WHEN {} >= {} AND {} < {} THEN {}\n'.format(var, array[i],var,array[i+1],math.sqrt(array[i]*array[i+1]))

    res+='ELSE NULL END '
    return res

def binLowEdgeFromArray(var, array,aliasName='binX'):
    res='CASE\n'
    for i in range(len(array)-1):
        res+='    WHEN {} >= {} AND {} < {} THEN {}\n'.format(var, array[i],var,array[i+1],array[i])
    res+='ELSE NULL END '
    return res

def binHighEdgeFromArray(var, array,aliasName='binX'):
    res='CASE\n'
    for i in range(len(array)-1):
        res+='    WHEN {} >= {} AND {} < {} THEN {}\n'.format(var, array[i],var,array[i+1],array[i+1])
    res+='ELSE NULL END '
    return res

def binIndex(nBins, firstBin, lastBin, var):
    # find in which bin is var
    binWidth=float(lastBin-firstBin)/nBins
    return " INTEGER(FLOOR( ( (" + var + ") - (" + str(firstBin) + "))/(" + str(binWidth) + ") )) "

def binCenter(nBins, firstBin, lastBin, var):
    # find in which bin is var
    binWidth=float(lastBin-firstBin)/nBins
    iBin=binIndex(nBins, firstBin, lastBin, var)
    return "({}+({}+0.5)*{})".format(firstBin,iBin,binWidth)

def binLowEdge(nBins, firstBin, lastBin, var):
    # find in which bin is var
    binWidth=float(lastBin-firstBin)/nBins
    iBin=binIndex(nBins, firstBin, lastBin, var)
    return "({}+({})*{})".format(firstBin,iBin,binWidth)

def binHighEdge(nBins, firstBin, lastBin, var):
    # find in which bin is var
    binWidth=float(lastBin-firstBin)/nBins
    iBin=binIndex(nBins, firstBin, lastBin, var)
    return "({}+({})*{}+1)".format(firstBin,iBin,binWidth)
    
def histCustomCommand( theCommand, requery=False):
    dirName=os.environ.get('HOME')+'/.bigQueryCached/'
    if requery is False:
        cachedResult=getHistDataFrame(dirName,theCommand)
        if cachedResult is not None:
            return cachedResult

    try:
        data = executeQuery(theCommand)
        df=pd.read_json(json.dumps(data))
        if not df.empty: 
            print 'saving'
            saveHistDataFrame(df,dirName,theCommand)

    except ValueError:
        print 'no json data'
        return None
    
    return df



def hist2DCustomCommand( nBinsX, firstBinX, lastBinX, nBinsY, firstBinY, lastBinY, theCommand):
    binWidth=float(lastBin-firstBin)/nBins

    try:
        data = executeQuery(theCommand)

        L=dict()
        L['Y']=dict()
        L['binX']=dict()
        L['binY']=dict()        

        L['Y']=[0 for i in range(nBins)]
        L['binX']=[firstBinX+0.5*binWidthX+i*binWidthX for i in range(nBinsX)]
        L['binY']=[firstBinY+0.5*binWidthY+i*binWidthY for i in range(nBinsY)]

        L = [[0 for x in range(nBinsX)] for y in range(nBinsY)] 

        for d in data:
            try:
                L['Y'][int(d['binX'])]=float(d['f0_'])
                L['Y'][int(d['binX'])]=float(d['f0_'])
            except:
                pass

        dataFrame=pd.DataFrame(L)
        return dataFrame
    except ValueError:
        print 'no json data'


class Hist:
    
    def __init__(self,df,nBinsX,firstBinX,lastBinX,nBinsY=None,firstBinY=None,lastBinY=None):
        self.df=df

        self.nBinsX=nBinsX
        self.firstBinX=firstBinX
        self.lastBinX=lastBinX
        self.nDimension=1

        if not nBinsY: return
            
        self.nBinsY=nBinsY
        self.firstBinY=firstBinY
        self.lastBinY=lastBinY        
        self.nDimension=2
        self.matrix = zeros([self.nBinsX,self.nBinsY])
        for index, row in self.df.iterrows():
            self.matrix[ row['binX'], row['binY'] ] = log10(row['f0_'])
        self.matrix=np.rot90(self.matrix)


    def sliceY(self,bin):
        df=self.df[self.df['binX']==bin]
        df=df.drop('binX', 1)
        df.columns=['binX','f0_']
        binWidthY=float(self.lastBinY-self.firstBinY)/self.nBinsY
        df['binX'] = self.firstBinY + (0.5+df['binX'])*binWidthY

        return Hist(df,self.nBinsY,self.firstBinY,self.lastBinY)

        
    def fit(self,fitFunction,plot=True,fittingRange=None,**kwargs):
        if self.df.empty:
            return
            
        X=self.df['binX'].tolist()
        Y=self.df['f0_'].tolist()

        if fittingRange is not None:
            index = (np.array(X) > fittingRange[0]) & (np.array(X) < fittingRange[1])
            X=np.array(X)[index]
            Y=np.array(Y)[index]

        popt, pcov = curve_fit(fitFunction,X,Y,**kwargs)

        if plot:
            print 'yolo'
            plt.figure()
            self.plot()
            plt.plot(X, fitFunction(X, *popt))

        return popt, pcov
                    

    def plot(self,**kwargs):
        if self.df.empty:
            return self
            
        if self.nDimension == 1:
            self.df.plot(x='binX',drawstyle='steps',**kwargs)
            plt.xticks(np.arange(self.firstBinX,self.lastBinX,float(self.lastBinX-self.firstBinX)/5))
            return self
        else:
            matshow(self.matrix,extent=[self.firstBinX,self.lastBinX,self.firstBinY,self.lastBinY],aspect="auto", **kwargs)
            plt.xticks(np.arange(self.firstBinX,self.lastBinX,float(self.lastBinX-self.firstBinX)/5))
            plt.yticks(np.arange(self.firstBinY,self.lastBinY,float(self.lastBinY-self.firstBinY)/5)) 
            return self

def hist( nBins, firstBin, lastBin, var, cut='', queryOption='', requery=False ):
    binWidth=float(lastBin-firstBin)/nBins
    theCommand="bq --format json query " + queryOption + " -n " + str(nBins+10) + " \"SELECT " + binCenter(nBins,firstBin,lastBin,var) + " as binX,COUNT(1) FROM " + bigQueryTable
    if cut :
        theCommand+=" WHERE (" + cut + ") "

    theCommand+=" GROUP BY binX"
    theCommand+=" HAVING binX >= " + str(firstBin) + " AND binX < (" + str(lastBin)+ " )"
    theCommand+=" ORDER BY binX\""

    df=histCustomCommand( theCommand, cacheName=cacheName )
    return Hist( df, nBins, firstBin, lastBin )
    

def hist2d( nBinsX, firstBinX, lastBinX, nBinsY, firstBinY, lastBinY, varX, varY, cut='', requery=False ):
    cacheName='hist2D_nBinsX={}_firstBinX={}_lastBinX={}_varX={}_nBinsY={}_firstBinY={}_lastBinY={}_varY={}_cut={}.pkl'.format(nBinsX,firstBinX,lastBinX,varX,nBinsY,firstBinY,lastBinY,varY,cut)


    binWidthX=float(lastBinX-firstBinX)/nBinsX
    binWidthY=float(lastBinY-firstBinY)/nBinsY

    theCommand="bq --format json query -n " + str(nBinsX*nBinsY+10) + " \"SELECT " + binIndex(nBinsX,firstBinX,lastBinX,varX) + " as binX, " + binIndex(nBinsY,firstBinY,lastBinY,varY) + " as binY,COUNT(1) FROM " + bigQueryTable
    if cut :
        theCommand+=" WHERE (" + cut + ")"

    theCommand+=" GROUP BY binX, binY"
    theCommand+=" HAVING binX >= 0 AND binX < (" + str(nBinsX) + ") AND binY >= 0 AND binY < (" + str(nBinsY) + ")"
    theCommand+=" ORDER BY binX,binY \""

    data = executeQuery(theCommand)
    df=pd.read_json(json.dumps(data))
    saveHistDataFrame(df,dirName,cacheName)

    
    histo=Hist(df,nBinsX, firstBinX, lastBinX, nBinsY, firstBinY, lastBinY)
    print 'the hist : {}'.format(histo)
    return histo
    # for d in data:
    #     L[int(d['binX'])*nBinsY+int(d['binY'])]=float(d['f0_'])
        
    # x=[firstBinX + (i/nBinsY)*binWidthX  for i in range(nBinsX*nBinsY)]
    # y=[firstBinY + (i%nBinsY)*binWidthY for i in range(nBinsX*nBinsY)]

    # plt.hist2d(x,y, bins=[nBinsX,nBinsY], range=[[firstBinX, lastBinX], [firstBinY, lastBinY]], weights=L)
    # plt.colorbar()

    # plt.show()
    # return hist


def setTable(newTableName):
    global bigQueryTable
    bigQueryTable=newTableName
    print 'global bigQueryTable set to : ' + bigQueryTable

def makeOldSelectionMask(cutList):
    data=executeQuery("bq --format json show " + bigQueryTable)
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
            print 'Cut ' + cut + 'not found !'
            return -1

    return theMask

def singleton(cls):
    obj = cls()
    # Always return the same object
    cls.__new__ = staticmethod(lambda cls: obj)
    # Disable __init__
    try:
        del cls.__init__
    except AttributeError:
        pass
        return cls

# @singleton
class SelStatusDescriptor:
    def __init__(self):
        data=executeQuery("bq --format json show " + bigQueryTable)
        for d in data['schema']['fields']:
            if d['name'] == 'selStatus':
                selStatusName=d['description']

        selStatusName=selStatusName.split(',')
        self.bitIndex=dict()
        for i in range(len(selStatusName)):
            self.bitIndex[selStatusName[i]]=i

def makeSelectionMask(cutList):
    selStatusDescriptor=SelStatusDescriptor()
    selMask=0 # has a 1 for every bit that has to be checked
    statusMask=0 # take the bit value for every bit that has to be checked
    for cut in cutList:
        val=1
        if cut[0]=='!':
            val=0
            cut=cut[1:]
            
        if cut in selStatusDescriptor.bitIndex:
            selMask+=1<<selStatusDescriptor.bitIndex[cut]
            statusMask+=val<<selStatusDescriptor.bitIndex[cut]
        else:
            print 'Cut ' + cut + ' not found !'
            return ''

    return ' ((selStatus^' + str(statusMask) + ')&' + str(selMask) + ')==0 '

def saveHistDataFrame(data,dirname,filename):
    try:
        f=open(dirname+'/'+str(filename.__hash__()),'wb')
        return pkl.dump(data,f)
    except IOError:
        try:
            os.mkdir(dirname)
        except:
            pass
    return None

def getHistDataFrame(dirname,filename):
    print 'loading : ' + filename
    try:
        f=open(dirname+'/'+str(filename.__hash__()))
        return pkl.load(f)
    except IOError:
        return None

