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

bigQueryTable="full_test.full_test"

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
        print rawOutput
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
    return dict()

# def ListDatasets(service, project):
#     try:
#         datasets = service.datasets()
#         list_reply = datasets.list(projectId=project).execute()
#         print 'Dataset list:'
#         pprint.pprint(list_reply)
#     except HTTPError as err:
#         print 'Error in ListDatasets:', pprint.pprint(err.content)
        
# Run a synchronous query, save the results to a table, overwriting the
# existing data, and print the first page of results.
# Default timeout is to wait until query finishes.
# def runSyncQuery (service, projectId, datasetId, timeout=0):
#     try:
#         print 'timeout:%d' % timeout
#         jobCollection = service.jobs()
#         queryData = {'query':'SELECT word,count(word) AS count FROM publicdata:samples.shakespeare GROUP BY word;','timeoutMs':timeout}

#     queryReply = jobCollection.query(projectId=projectId, body=queryData).execute()

#     jobReference=queryReply['jobReference']

#     # Timeout exceeded: keep polling until the job is complete.
#     while(not queryReply['jobComplete']):
#         print 'Job not yet complete...'
#         queryReply = jobCollection.getQueryResults(
#             projectId=jobReference['projectId'],
#             jobId=jobReference['jobId'],
#             timeoutMs=timeout).execute()

#     # If the result has rows, print the rows in the reply.
#     if('rows' in queryReply):
#         print 'has a rows attribute'
#         printTableData(queryReply, 0)
#         currentRow = len(queryReply['rows'])

#       # Loop through each page of data
#       while('rows' in queryReply and currentRow < queryReply['totalRows']):
#           queryReply = jobCollection.getQueryResults(
#               projectId=jobReference['projectId'],
#               jobId=jobReference['jobId'],
#               startIndex=currentRow).execute()
#           if('rows' in queryReply):
#               printTableData(queryReply, currentRow)
#               currentRow += len(queryReply['rows'])

#   except AccessTokenRefreshError:
#       print ("The credentials have been revoked or expired, please re-run"
#              "the application to re-authorize")

#   except HttpError as err:
#       print 'Error in runSyncQuery:', pprint.pprint(err.content)

#   except Exception as err:
#       print 'Undefined error' % err                                                                                                                  
      
def executeQueryNew(theCommand):
    client = bq.Client.Get()
    tableid = client.Query(theCommand)['configuration']['query']['destinationTable']
    table = client.ReadTableRows(tableid)
    print table
    return table

def bin(nBins, firstBin, lastBin, var):
    # find in which bin is var
    binWidth=float(lastBin-firstBin)/nBins
    return " INTEGER(FLOOR( ( (" + var + ") - (" + str(firstBin) + "))/(" + str(binWidth) + ") )) "

def histCustomCommandOld( nBins, firstBin, lastBin, theCommand):
    binWidth=float(lastBin-firstBin)/nBins

    try:
        data = executeQuery(theCommand)

        L=[0 for i in range(nBins)]
        
        for d in data:
            try:
                L[int(d['binX'])]=float(d['f0_'])
            except:
                pass
            
        x=[firstBin + i*binWidth for i in range(nBins)]
        hist, bins = np.histogram(x, range=(firstBin,lastBin),bins=nBins, weights=L)

        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        #    plt.yscale('log')
        plt.bar(center, hist, width=width)
        plt.grid()
        return hist, center, width
    except ValueError:
        print 'no json data'

def histCustomCommand( nBins, firstBin, lastBin, theCommand):
    binWidth=float(lastBin-firstBin)/nBins

    try:
        data = executeQuery(theCommand)

        L=[0 for i in range(nBins)]

        for d in data:
            try:
                L[int(d['binX'])]=float(d['f0_'])
            except:
                pass

        dataFrame=pd.DataFrame(np.array(L),columns=['Y'])
        return dataFrame
    except ValueError:
        print 'no json data'

def hist( nBins, firstBin, lastBin, var, cut='', queryOption='' ):
    binWidth=float(lastBin-firstBin)/nBins
    theCommand="bq --format json query " + queryOption + " -n " + str(nBins+10) + " \'SELECT " + bin(nBins,firstBin,lastBin,var) + " as binX,COUNT(1) FROM " + bigQueryTable
    if cut :
        theCommand+=" WHERE (" + cut + ") "

    theCommand+=" GROUP BY binX"
    theCommand+=" HAVING binX >= 0 AND binX < (" + str(nBins)+ " )"
    theCommand+=" ORDER BY binX\'"

    # h = histQueryFactory.HistQueryFactory()
    # h.add_variable(var,nBins,firstBin,lastBin)
    # h.add_condition(cut)
    #theSelectClause=str(h)
    #theCommand="bq --format json query " + queryOption + " -n " + str(nBins+10) + " \'" + theSelectClause + " \'" 
    
    return histCustomCommand( nBins, firstBin, lastBin, theCommand )

def hist2d( nBinsX, firstBinX, lastBinX, nBinsY, firstBinY, lastBinY, varX, varY, cut='' ):
    binWidthX=float(lastBinX-firstBinX)/nBinsX
    binWidthY=float(lastBinY-firstBinY)/nBinsY

    theCommand="bq --format json query -n " + str(nBinsX*nBinsY+10) + " \'SELECT " + bin(nBinsX,firstBinX,lastBinX,varX) + " as binX, " + bin(nBinsY,firstBinY,lastBinY,varY) + " as binY,COUNT(1) FROM " + bigQueryTable
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
    return hist


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

def makeSelectionMask(cutList):
    data=executeQuery("bq --format json show " + bigQueryTable)
    for d in data['schema']['fields']:
        if d['name'] == 'selStatus':
            selStatusName=d['description']
    
    selStatusName=selStatusName.split(',')
    bitIndex=dict()
    for i in range(len(selStatusName)):
        bitIndex[selStatusName[i]]=i

    selMask=0 # has a 1 for every bit that has to be checked
    statusMask=0 # take the bit value for every bit that has to be checked
    for cut in cutList:
        val=1
        if cut[0]=='!':
            val=0
            cut=cut[1:]
            
        if cut in bitIndex:
            selMask+=1<<bitIndex[cut]
            statusMask+=val<<bitIndex[cut]
        else:
            print 'Cut ' + cut + ' not found !'
            return ''

    return ' ((selStatus^' + str(statusMask) + ')&' + str(selMask) + ')==0 '

def saveHist(center,hist,width,filename):
    data=[hist,center,width]
    f=open(filename,'wb')
    return pkl.dump(data,f)

def getHist(filename):
    f=open(filename)
    data=pkl.load(f)
    return plt.bar(data[0], data[1], width=data[2])
