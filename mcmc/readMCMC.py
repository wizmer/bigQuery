import MCMC as mcmc
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from pylab import hist, show
import gzip
import time
import sys
import getopt

class Reader:
    verbose = False
    
    def __init__(self,filename,dirname='./',maxChunk=0):
        print 'maxchunk : '+str(maxChunk)
        self.filename=filename
        self.dirname=dirname
        self.theFileName=dirname+'/thread0_'+filename
        self.maxChunk=maxChunk
        f=open(self.theFileName,'r')
        self.metadata=pickle.load(f)
        f.close()
        print self.metadata

        self.burnInLength=0
        self.correlationLength=1

        self.totalCount=np.zeros((self.metadata['nVar'],100))
        self.heatmap=np.zeros((self.metadata['nVar'],self.metadata['nVar'],100,100))
        self.extent=np.zeros((self.metadata['nVar'],self.metadata['nVar'],4))

        self.bins=np.zeros((self.metadata['nVar'],101))
        self.firstHisto=True
        self.chunkNumber=0

        self.t0=time.clock()
        # self.array=[0]*self.metadata['nVar']

        self.load()
        self.plot()

    def loadAll(self):
        self.load1D()
        self.load2D()

    def load1D(self):
        if self.verbose: print '1d'

        if self.firstHisto:
            if self.verbose: print 'first'
            for iVar in range(self.metadata['nVar']):
                print '1D var #'+str(iVar)
                counts, self.bins[iVar], totalPatch=hist(self.df[iVar],bins=100)
                self.totalCount[iVar]+=counts
        else:
            if self.verbose: print 'not first'
            for iVar in range(self.metadata['nVar']):
                counts, b, totalPatch=hist(self.df[iVar],bins=self.bins[iVar])
                self.totalCount[iVar]+=counts
                if iVar == 0:
                    if self.verbose: print counts

    def load1DOld(self):
        if self.firstHisto:
            if self.verbose: print 'first'
            for iVar in range(self.metadata['nVar']):
                if self.verbose: print '1D var #'+str(iVar)
                counts, self.bins[iVar], totalPatch=hist(self.array[iVar],bins=100)
                self.totalCount[iVar]+=counts
        else:
            if self.verbose: print 'not first'
            for iVar in range(self.metadata['nVar']):
                counts, b, totalPatch=hist(self.array[iVar],bins=self.bins[iVar])
                self.totalCount[iVar]+=counts
                if iVar == 0:
                    if self.verbose: print counts

    def load2D(self):
        if self.verbose: print '2d'
        if self.verbose: print self.firstHisto
        if self.firstHisto:
            for iVar in range(self.metadata['nVar']):
                for jVar in range(self.metadata['nVar']):
                    if self.verbose: print '2d plots: {},{}'.format(iVar,jVar)
                    self.heatmap[iVar][jVar], xedges, yedges = np.histogram2d(self.df[iVar],self.df[jVar],bins=100)
                    self.extent[iVar][jVar] = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        else:
            size=len(self.df.index)
            if self.verbose: print 'size : {}'.format(self.metadata['nVar'])
            for iVar in range(self.metadata['nVar']):
                for jVar in range(self.metadata['nVar']):
                    if self.verbose: print '2d plots: {},{}, nentries: {}'.format(iVar,jVar,size)
                    h, xedges, yedges = np.histogram2d(self.df[iVar],self.df[jVar],bins=100,range=[[self.extent[iVar][jVar][0],self.extent[iVar][jVar][1]],[self.extent[iVar][jVar][2],self.extent[iVar][jVar][3]]])
                    self.heatmap[iVar][jVar]+=h


    def load2DOld(self):
        if self.firstHisto:
            for iVar in range(self.metadata['nVar']):
                for jVar in range(self.metadata['nVar']):
                    print '2d plots: {},{}'.format(iVar,jVar)
                    self.heatmap[iVar][jVar], xedges, yedges = np.histogram2d(self.array[iVar],self.array[jVar],bins=100)
                    self.extent[iVar][jVar] = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            else:
                for iVar in range(self.metadata['nVar']):
                    for jVar in range(self.metadata['nVar']):
                        print '2d plots: {},{}, nentries: {}'.format(iVar,jVar,len(self.array[iVar]))
                        h, xedges, yedges = np.histogram2d(self.array[iVar],self.array[jVar],bins=100,range=[[self.extent[iVar][jVar][0],self.extent[iVar][jVar][1]],[self.extent[iVar][jVar][2],self.extent[iVar][jVar][3]]])
                        self.heatmap[iVar][jVar]+=h

    def loadOld(self):
        for iFile in range(4):
            self.theFileName=self.dirname+'/thread{}_'.format(iFile)+self.filename
            print 'opening '+self.theFileName
            try:
                f=open(self.theFileName,'r')
            except:
                print "Error loading chunk:", sys.exc_info()[0]
                break
            
            test=pickle.load(f) #skip metadata
            remainingBurnInLength=self.burnInLength
            while True:
                try:
                    print 'loading chunkNumber : '+str(self.chunkNumber)
                    if self.maxChunk > 0 and self.chunkNumber >= self.maxChunk:
                        break
                    self.chunk=pickle.load(f)
                    self.chunkNumber+=1
                    print 'chunk loaded'

                    # skip chunk before burn-in length
                    if remainingBurnInLength >= self.chunk['numberOfSteps']:
                        print 'skip chunk'
                        remainingBurnInLength-=self.chunk['numberOfSteps']
                        if remainingBurnInLength < 0: remainingBurnInLength=0
                        continue
                except:
                    print "Error loading chunk:", sys.exc_info()[0]
                    break

                print 'start array making'
                for iVar in range(self.metadata['nVar']):
                    self.array[iVar]=np.array(self.chunk['trace'])[:,iVar][remainingBurnInLength:self.chunk['numberOfSteps']:self.correlationLength]
                print 'end array making'

                try:
                    self.loadAll()
                    self.firstHisto=False
                except:
                    print "histo creation failed:", sys.exc_info()[0]
                    pass

                remainingBurnInLength-=self.chunk['numberOfSteps']
                if remainingBurnInLength < 0: remainingBurnInLength=0

                del self.chunk
            f.close()


    def load(self):
        for iFile in range(4):
            self.theFileName=self.dirname+'/thread{}_'.format(iFile)+self.filename
            print 'opening '+self.theFileName
            try:
                f=open(self.theFileName,'r')
            except EOFError:
                break
            except:
                print "Error reading file:", sys.exc_info()[0]
                break
            
            test=pickle.load(f) #skip metadata
            remainingBurnInLength=self.burnInLength
            while True:
                try:
                    print 'loading chunkNumber : '+str(self.chunkNumber)
                    if self.maxChunk > 0 and self.chunkNumber >= self.maxChunk:
                        break
                    self.df=pickle.load(f)
                    self.chunkNumber+=1
                    print 'chunk loaded'

                    # skip chunk before burn-in length
                    # if remainingBurnInLength >= self.chunk['numberOfSteps']:
                    #     print 'skip chunk'
                    #     remainingBurnInLength-=self.chunk['numberOfSteps']
                    #     if remainingBurnInLength < 0: remainingBurnInLength=0
                    #     continue
                except EOFError:
                    break
                except:
                    print "Error loading chunk:", sys.exc_info()[0]
                    break

                try:
                    self.loadAll()
                    self.firstHisto=False
                except:
                    print "histo creation failed:", sys.exc_info()[0]
                    pass

                # remainingBurnInLength-=self.chunk['numberOfSteps']
                # if remainingBurnInLength < 0: remainingBurnInLength=0
            f.close()
            
    def plot(self):
        print time.clock()-self.t0
        print 'plotting'
        plt.figure('Correlation plots',figsize=(15,20))
        for iVar in range(self.metadata['nVar']):
            print 'iVar: '+str(iVar)
            for jVar in range(self.metadata['nVar']):
                print 'jVar: '+str(jVar)
                plotNum=iVar*self.metadata['nVar']+jVar+1
                plt.subplot(self.metadata['nVar'],self.metadata['nVar'],plotNum)
                if iVar==jVar:
                    plt.plot(self.bins[iVar][:-1],self.totalCount[iVar],drawstyle='steps')
                    if len(self.metadata['realValues']) > 0:
                        plt.axvline(self.metadata['realValues'][iVar])
                elif iVar>jVar:
                    plt.imshow(self.heatmap[jVar][iVar], extent=self.extent[jVar][iVar],aspect="auto")
                    plt.title('titre {},{}'.format(iVar,jVar))

        plt.savefig(self.filename+'.png')
        plt.savefig(self.filename+'.eps')



class histo1D:
    def __init__(self,metadata,bins,counts):
        self.metadata=metadata
        self.bins=bins
        self.counts=counts
        
class Reader1D(Reader):
    def loadAll(self):
        print 'loading all'
        self.load1D()

    def plot(self):
        print time.clock()-self.t0
        print 'plotting'
        plt.figure('Correlation plots')
        f=open(self.filename+'.pkl','w')
        pickle.dump(self.metadata,f)
        histos=[]
        for iVar in range(self.metadata['nVar']):
            print 'iVar: '+str(iVar)
            plt.subplot(self.metadata['nVar'],1,iVar+1)
            h=histo1D(self.metadata,self.bins[iVar],self.totalCount[iVar])
            histos.append(h)
            plt.plot(self.bins[iVar][:-1],self.totalCount[iVar],drawstyle='steps')
            if len(self.metadata['realValues']) > 0:
                plt.axvline(self.metadata['realValues'][iVar])
                
        pickle.dump(histos,f,2)

        plt.savefig(self.filename+'.png')
        plt.savefig(self.filename+'.eps')

def main(argv):
    filename='test.pkl'
    dirname='./'
    theMaxChunk=0
    print argv
    try:
        opts, args = getopt.getopt(argv,"f:d:c:")
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    
    for opt, arg in opts:
        print 'option {},{}'.format(opt,arg)
        if opt == '-f':
            filename=arg
        elif opt == '-d':
            dirname=arg
        elif opt == '-c':
            theMaxChunk=int(arg)

    print filename
    Reader1D(filename,dirname,maxChunk=theMaxChunk)
    show()

if __name__ == "__main__":
    main(sys.argv[1:])


