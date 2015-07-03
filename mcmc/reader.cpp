#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TKey.h"

#include "generalUtils.hpp"
#include "Stack.hpp"

constexpr bool monitorTime = false;

enum ReadMode{
    kFull,
    kSparse
};

class ReaderBase{
public:
    ReaderBase(std::string _filename){
        filename = _filename;
        maxSteps = 0;
        nPlotPerCanvas = 4;
        isToyMC = false;
        
        isAlreadyAnalysed = generalUtils::fileExists(filename+"/.ana");
        init();

    }

    void setMaxSteps(int _maxSteps){
        maxSteps = _maxSteps;
    }

    void setCorrelationLength(int _correlationLength){
        correlationLength = _correlationLength;
    }


    void rePlotOffset(){
        std::vector<std::string> files = generalUtils::splitIntoLines( generalUtils::exec("find "+ filename + " | grep fullCan"));
        for(int i = 0;i<files.size();i++){
            std::cout << "files[i] : " << files[i] << std::endl;
            TFile* file = new TFile( files[i].c_str() );
            TKey*k = (TKey*)(file -> GetListOfKeys() -> At(0));
            std::cout << "k : " << k << std::endl;
            TCanvas* can = (TCanvas*) k -> ReadObj();
            can -> Draw();
	}
    }

    void setMaxChunk(int _maxChunk ){ maxChunk = _maxChunk; }
    void setFirstChunk(int _firstChunk ){ firstChunk = _firstChunk; }

    void rePlotFull(){

    }

    void rePlotAll(){
        rePlotOffset();
        rePlotFull();
    }

    void readAll(){
        if( isAlreadyAnalysed ) return rePlotAll();
        // Extract number of variables

        int nChunks = generalUtils::stringTo<int>( generalUtils::exec("ls -lrt " + filename + " | grep -o -E \"chunk[0-9]+\" | cut -c 6- | sort | uniq | wc -l") );
        std::cout << "nChunks : " << nChunks << std::endl;
        int lastChunk = nChunks;
        if( maxChunk > -1 && maxChunk < nChunks - firstChunk) lastChunk = maxChunk + firstChunk;

        prepareHisto(firstChunk,lastChunk);

        ReadMode readMode;
        readMode = ReadMode::kFull;
        //readMode = ReadMode::kSparse;
        
        std::cout << "lastChunk : " << lastChunk << std::endl;
        std::cout << "maxChunk : " << maxChunk << std::endl;
        std::cout << "nChunks : " << nChunks << std::endl;
        std::cout << "firstChunk : " << firstChunk << std::endl;
        for(int jChunk = firstChunk; jChunk<lastChunk; jChunk++){
            for(int iVar = 0; iVar<nVar+1 ;iVar++){
                //std::cout << "jChunk : " << jChunk << "\tiVar : " << iVar << std::endl;
                unsigned int nEntries = readBinary( getChunkName(jChunk,iVar),  iVar, readMode );
                fillHistos(iVar, nEntries);
            }
                
        }
        plot();
    }
    
    void prepareHisto(int firstChunk, int lastChunk){
        std::cout << "prepareHisto" << std::endl;
        float val = 0;
        for(int iVar = 0; iVar<nVar+1 ;iVar++){
            float sum = 0;
            float sum2 = 0;
            int n;
            int count = 0;
            for(int jChunk = firstChunk; jChunk<lastChunk; jChunk++){
                // Rough estimation of the mean and variance of the variable to create the histogram
                unsigned int nEntries = 0;
                std::ifstream f( getChunkName(jChunk,iVar).c_str(), std::ios::binary );
                f.read((char*)&nEntries, sizeof(nEntries));

                n = 100;
                if( nEntries < n ) n = nEntries;

                for(int i = 0; i<n;i++){
                    if( maxSteps > 0 && nEntries > maxSteps ) nEntries = maxSteps;
                    f.seekg(sizeof(nEntries)+i*nEntries/n*sizeof(float));
                    f.read((char*)&val, sizeof(float));
                    sum +=  val;
                    sum2 += pow(val, 2);
                    count++;
                }
            }

            float mean = sum/(float)count;
            float sigma = sqrt( sum2/(float)count - mean*mean);
            std::cout << "sigma : " << sigma << std::endl;
            borne[iVar].first  = mean - 10*sigma;
            borne[iVar].second = mean + 10*sigma;
            
            initHisto(iVar);

        }

        std::cout << "End prepareHisto" << std::endl;
    }

protected:
    std::string filename;
    virtual void plot() = 0;

    virtual void fillHistos(int, unsigned int) = 0;
        
    void init(){
        readMetaData();
        trace = new float*[(int)metadata["nVar"]+1]; // Last index is Likelihood value
        firstRead = new bool[(int)metadata["nVar"]+1];
        borne = new std::pair<int,int>[nVar+1];
        for(int i = 0;i<metadata["nVar"]+1;i++){
            trace[i] = new float[(int)metadata["chunkSize"]];
            firstRead[i] = true;
            borne[i] = std::pair<int,int>();
        }
        firstChunk = 0;
        maxChunk = -1;
        correlationLength = 1;
        nBins = 500;
    }
    
    std::string getChunkName(int iChunk, int iVar){
        if( iVar == nVar ) return filename+"/lik_chunk" + generalUtils::toString(iChunk) + ".bin";
        return filename + "/par" + generalUtils::toString(iVar) + "_chunk" + generalUtils::toString(iChunk) + ".bin";
    }

    
    std::map< std::string, float>  metadata;
    void readMetaData(){
        std::ifstream f( (filename+"/metadata.txt").c_str());
        while( ! f.eof() ){
            std::string key;
            float value;
            f >> key >> value;
            if( key == "" ) break;
            metadata[key] = value;
            std::cout << "key : " << key << std::endl;
            std::cout << "value : " << value << std::endl;

        }
        nVar = metadata["nVar"];
        isToyMC = metadata["isToyMC"];
        if( maxSteps > 0 && metadata["nStep"] > maxSteps ) metadata["nStep"] = maxSteps;
    }
    
    virtual void initHisto(int iParam){
        std::cout << iParam << "\t" << borne[iParam].first << "\t" << borne[iParam].second << std::endl;
        h[iParam] = new TH1F(Form("h_%i",iParam),Form("h_%i",iParam), nBins,borne[iParam].first,borne[iParam].second);
        
        for(int j = 0;j<iParam;j++){
            h2[ std::pair<int,int>(iParam,j) ] = new TH2F(Form("h2_%i_%i",iParam,j),Form("h2_%i_%i",iParam,j), 100,borne[iParam].first, borne[iParam].second, 100, borne[j].first, borne[j].second);
        }
    }

    unsigned int readFullBinary(std::string fname, int iParam ){
        std::clock_t start;
        if(monitorTime) start = std::clock();

        std::ifstream f(fname.c_str(), std::ios::binary );

        float sum = 0;

        unsigned int nEntries;

        f.read((char*)&nEntries, sizeof(unsigned int));
        if( maxSteps > 0 && nEntries > maxSteps ) nEntries = maxSteps;
        f.read((char*)trace[iParam], sizeof(float)*nEntries);
        if(monitorTime) std::cout << "Time readFullBinary : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        return nEntries;
    }

    unsigned int readSparseBinary(std::string fname, int iParam ){
        std::clock_t start;
        if(monitorTime) start = std::clock();

        std::ifstream f(fname.c_str(), std::ios::binary );

        float sum = 0;

        unsigned int nEntries;

        f.read((char*)&nEntries, sizeof(unsigned int));
        //        if( maxSteps > 0 && nEntries > maxSteps ) nEntries = maxSteps;
        
        for(int i = 0;i<nEntries/(float)correlationLength;i++){
            f.seekg(sizeof(nEntries)+i*correlationLength*sizeof(float));
            f.read((char*)&trace[iParam][i*correlationLength], sizeof(float));
	}

        if(monitorTime) std::cout << "Time readSparseBinary : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        return nEntries/(float)correlationLength;
    }

    unsigned int readBinary(std::string fname, int iParam, ReadMode readMode ){
        if( readMode == kFull ) return readFullBinary( fname, iParam);
        if( readMode == kSparse ) return readSparseBinary( fname, iParam);
    }

    virtual void smartPlot(const std::vector<TH1*> &histo){
        static int smartPlotCounter = 0;
        int nHistos = histo.size();
        if( nPlotPerCanvas > nHistos ) nPlotPerCanvas = nHistos;
        for( int iCan = 0; iCan < nHistos/nPlotPerCanvas+1*(nHistos%nPlotPerCanvas != 0); iCan++){
            TCanvas* can = new TCanvas(Form("can%i_%i",smartPlotCounter,iCan),Form("can%i_%i",smartPlotCounter,iCan),800,800);
            can -> Divide(1, nPlotPerCanvas );
            for(int iPlot = 0;iPlot < nPlotPerCanvas;iPlot++){
                int iVar = iPlot+iCan*nPlotPerCanvas;
                if( iVar >= nHistos ) return;
                can -> cd( 1+iPlot );
                Stack* st = new Stack();
                st -> push_back( histo[iVar] );
                st -> draw(gPad);
                // histo[iVar] -> Draw("colz");
            }
        }
        smartPlotCounter++;
    }

    float **trace;
    std::map< int, TH1* > h;
    std::map< std::pair<int,int>, TH2* > h2;
    int nVar;
    int maxSteps;
    int nPlotPerCanvas;
    bool *firstRead;
    std::pair<int,int> *borne;
    bool isToyMC;
    bool isAlreadyAnalysed;
    int firstChunk;
    int maxChunk;
    int correlationLength;
    int nBins;
};

class FullReader: public ReaderBase{
public:
    FullReader(std::string _filename): ReaderBase(_filename){};

protected:
    virtual void plot(){
        std::string resultDir = filename+"Ana/";
        generalUtils::makeFolder(resultDir);
        generalUtils::exec("touch " + resultDir + ".ana");

        std::map< int, Stack* > stackMap;
        TGraphErrors* grOffset[2];

        enum Species{
            P,
            D
        } species;

        int graphCounter[2] = {0,0};

        if(isToyMC){
            grOffset[P] = new TGraphErrors();
            grOffset[D] = new TGraphErrors();
        }

        if( nPlotPerCanvas > nVar ) nPlotPerCanvas = nVar;
        for( int iCan = 0; iCan < nVar/nPlotPerCanvas+1*(nVar%nPlotPerCanvas != 0); iCan++){
            TCanvas* can = new TCanvas(Form("fullCan_%i",iCan));
            can -> Divide( nPlotPerCanvas, nPlotPerCanvas );
            for(int iPlot = 0;iPlot < nPlotPerCanvas;iPlot++){
                int iVar = iPlot+iCan*nPlotPerCanvas;
                if( iVar >= nVar ) break;
                for(int jPlot = 0; jPlot <= iPlot; jPlot++){
                    int jVar = jPlot+iCan*nPlotPerCanvas;
                    can -> cd( 1+iPlot*nPlotPerCanvas+jPlot );
                    
                    stackMap[iVar] = new Stack(Form("Param #%i",iVar));
                    if( iVar==jVar ){
                        stackMap[iVar] -> push_back( h[iVar] );
                        std::cout << "h[iVar] : " << h[iVar]->GetName() << std::endl;
                        
                        float mean = h[iVar] -> GetMean();
                        float sigma = h[iVar] -> GetStdDev();
                        TF1* f = new TF1("f","[0]*exp(-pow(([1]-x)/[2],2))", mean - 2*sigma, mean + 2*sigma);
                        f -> SetParameter(0, h[iVar] -> GetMaximum() );
                        f -> SetParameter(1, mean );
                        f -> SetParameter(2, sigma );
                        h[iVar] -> Fit(f,"R0");
                        stackMap[iVar] -> push_back(f);
                        if(isToyMC){
                            stackMap[iVar] -> pushVerticalLine( metadata[ Form("realValues_%i",iVar) ] );
                            //if( std::abs( f-> GetParameter(1) - metadata[ Form("realValues_%i",iVar) ] ) < 1000 ){
                            if( true ){
                                Species s = iVar<nVar/2?P:D;
                                grOffset[s] -> SetPoint(graphCounter[s],iVar-(nVar/2)*(iVar>=nVar/2), f-> GetParameter(1) - metadata[ Form("realValues_%i",iVar) ]);
                                grOffset[s] -> SetPointError(graphCounter[s]++,0, f-> GetParameter(2));
                            }
                        }
                        
                    } else stackMap[iVar] -> push_back( h2[ std::pair<int, int>(iVar,jVar) ] );
                    stackMap[iVar] -> draw(gPad);
                }
            }
            can -> SaveAs( (resultDir+can->GetName()+".root").c_str() );
        }

        std::cout << "isToyMC : " << isToyMC << std::endl;
        if( isToyMC ){
            Stack* st = new Stack("offset");
            st -> push_back(grOffset[P],"proton");
            st -> push_back(grOffset[D],"deuton");
            st -> draw();
            st -> write((resultDir+"offset"),"root");
        }
    }

    virtual void fillHistos(int iParam, unsigned int nEntries){
        std::clock_t start;
        if( monitorTime ) start = std::clock();

        for(unsigned int j = 0;j < nEntries/correlationLength ;j++){
            h[iParam] -> Fill( trace[iParam][j*correlationLength] );
            //std::cout << "trace["<<iParam<<"]["<< j << "] = " << trace[iParam][j] << std::endl;
        }

        for(int j = iParam-3;j<iParam;j++){
            if( j < 0 ) continue;
            for(int k = 0;k < nEntries/correlationLength ;k++) h2[ std::pair<int,int>(iParam,j) ] -> Fill( trace[iParam][k*correlationLength], trace[j][k*correlationLength] );
        }

        if( monitorTime ) std::cout << "Time FillHisto : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    }
    

};

class LikelihoodTracer: public ReaderBase{
public:
    LikelihoodTracer(std::string _filename): ReaderBase(_filename){};

protected:
    virtual void plot(){
        // TCanvas* can = new TCanvas();
        // Stack* st = new Stack(Form("Param #%i",0));
        // st -> push_back( h[metadata["nVar"]] );
        // st -> draw(gPad);

        std::vector<TH1*> likelihoodVsParam;
        for(int iParam = 0;iParam<nVar;iParam++){
            likelihoodVsParam.push_back(h2[ std::pair<int,int>(nVar,iParam) ]);
        }
        smartPlot( likelihoodVsParam );

    }

    virtual void fillHistos(int iParam, unsigned int nEntries){
        if(iParam != nVar ) return;

        for(int j = 0;j<nVar;j++){
            for(unsigned int k = 0;k < nEntries ;k++) h2[ std::pair<int,int>(iParam,j) ] -> Fill( trace[iParam][k], trace[j][k] );
        }

        for(unsigned int k = 0;k < nEntries ;k++) h[ nVar ] -> Fill( trace[nVar][k] );
    }

    void initHisto(int iParam){
        if(iParam < nVar) return;
        h[iParam] = new TH1F(Form("h_%i",iParam),Form("h_%i",iParam), nBins,borne[iParam].first,borne[iParam].second);
        
        for(int j = 0;j<nVar;j++){
            if( firstRead[j] == false ){
                h2[ std::pair<int,int>(iParam,j) ] = new TH2F(Form("h2_%i_%i",iParam,j),Form("h2_%i_%i",iParam,j), 100,borne[iParam].first, borne[iParam].second, 100, borne[j].first, borne[j].second);
            }
        }
    }
};

class Tracer: public ReaderBase{
public:
    Tracer(std::string _filename): ReaderBase(_filename){};

protected:
    virtual void plot(){
        std::vector<TH1*> traces;
        for(int iParam = 0;iParam<nVar;iParam++){
            traces.push_back(h[iParam]);
        }

        smartPlot( traces );

    }

    virtual void fillHistos(int iParam, unsigned int nEntries){
        if(iParam == nVar ) return;

        for(unsigned int k = 0;k < nEntries ;k++) h[ iParam ] -> Fill( k, trace[iParam][k] );
        //std::cout << "trace[iParam][metadata[nStep]-1] : " << trace[iParam][(int)(metadata["nStep"]-1)] << std::endl;
    }

    void initHisto(int iParam){
        if(iParam == nVar) return;
        h[iParam] = new TH2F(Form("h_%i",iParam),Form("h_%i",iParam), 1000,0,metadata["nStep"],100,borne[iParam].first,borne[iParam].second);
    }
};

class LastValue: public ReaderBase{
public:
    LastValue(std::string _filename): ReaderBase(_filename){};

protected:
    std::vector<float> initialCondition;
    virtual void plot(){
        for(int i = 0;i<initialCondition.size()/2;i++) std::cout << initialCondition[i] << "\t" << initialCondition[i+initialCondition.size()/2] << std::endl;
    }

    virtual void fillHistos(int iParam, unsigned int nEntries){
        initialCondition.push_back(trace[iParam][nEntries-1]);
    }

    void initHisto(int iParam){
    }
};


int main(int argc, char** argv){
    TApplication app("app",&argc,argv,0,-1);

    int c;
    std::string dirname = "test";
    std::cout << "argv[1] : " << argv[1] << std::endl;
    int maxChunk = -1, firstChunk = 0, correlationLength = 1;
    if( argc > 1 ) dirname = argv[1];

    ReaderBase *r = NULL;

    while((c =  getopt(argc, argv, "t:c:f:l:")) != EOF)
        {
            switch (c)
                {
                case 't':
                    {
                        std::string a(optarg);
                        if( a == "full" ){
                            r = new FullReader(dirname);
                        } else if( a == "lik" ){
                            r = new LikelihoodTracer(dirname);
                        } else if( a == "trace" ){
                            r = new Tracer(dirname);
                        } else if( a == "value" ){
                            r = new LastValue(dirname);
                        }
                        break;
                    }
                case 'c':
                    {
                        maxChunk = generalUtils::stringTo<int>(optarg);
                        break;
                    }
                case 'f':
                    {
                        firstChunk = generalUtils::stringTo<int>(optarg);
                        break;
                    }
                case 'l':
                    {
                        correlationLength = generalUtils::stringTo<int>(optarg);
                        break;
                    }
                }
        }
    
    std::clock_t start = std::clock();

    if( r == NULL ) r = new FullReader(dirname);

    //    r.setMaxSteps(1000);
    r -> setMaxChunk(maxChunk);
    r -> setFirstChunk(firstChunk);
    r -> setCorrelationLength(correlationLength);
    r -> readAll();


    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    app.Run();
    return 0;
}
