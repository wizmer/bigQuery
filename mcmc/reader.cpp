#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "generalUtils.hpp"
#include "Stack.hpp"

class ReaderBase{
  std::string filename;
public:
  ReaderBase(std::string _filename){
    filename = _filename;
    maxSteps = 0;
    nPlotPerCanvas = 4;
  }

  void setMaxSteps(int _maxSteps){
    maxSteps = _maxSteps;
  }

  void readAll(){
    init();
    // Extract number of variables

    int nChunks = generalUtils::stringTo<int>( generalUtils::exec("ls -lrt " + filename + " | grep -o \"chunk[0-9]*\" | cut -c 6- | sort | uniq | wc -l") );
    std::cout << "nChunks : " << nChunks << std::endl;

    for(int jChunk = 0; jChunk<nChunks; jChunk++){
      for(int iVar = 0; iVar<nVar+1 ;iVar++){
        std::string chunkName = Form("%s/par%i_chunk%i.bin", filename.c_str(), iVar, jChunk);
        if( iVar == nVar ) chunkName = Form("%s/lik_chunk%i.bin", filename.c_str(), jChunk);
        readBinary( chunkName, iVar );
        fillHistos(chunkName,iVar);
      }
                
    }

    plot();
  }
    
protected:
  virtual void plot() = 0;

  virtual void fillHistos(std::string, int) = 0;
        
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

    if( maxSteps > 0 && metadata["nStep"] > maxSteps ) metadata["nStep"] = maxSteps;
  }
    
  virtual void initHisto(int iParam){
    h[iParam] = new TH1F(Form("h_%i",iParam),Form("h_%i",iParam), 100,borne[iParam].first,borne[iParam].second);
        
    for(int j = 0;j<nVar;j++){
      if( firstRead[j] == false ){
        h2[ std::pair<int,int>(iParam,j) ] = new TH2F(Form("h2_%i_%i",iParam,j),Form("h2_%i_%i",iParam,j), 100,borne[iParam].first, borne[iParam].second, 100, borne[j].first, borne[j].second);
      }
    }
  }

  void readBinary(std::string fname, int iParam){
    std::ifstream f(fname.c_str(), std::ios::binary );

    float sum = 0;

    unsigned int N;

    f.read((char*)&N, sizeof(unsigned int));
    if( maxSteps > 0 && N > maxSteps ) N = maxSteps;
    f.read((char*)trace[iParam], sizeof(float)*N);

    // Rough estimation of the mean and variance of the variable to create the histogram
    if( firstRead[iParam] ){
      int n = 100;
      if( N < n ) n = N;
      float sum = 0;
      float sum2 = 0;
      for(int i = 0; i<n;i++){
        int index = (int)(i/(float)n * N);
        sum +=  trace[iParam][index];
        //        std::cout << trace[iParam][index] << "\t" << iParam << "\t" << index << std::endl;
        sum2 += pow(trace[iParam][index], 2);
        //std::cout << "trace["<<iParam<<"]["<< index << "] = " << trace[iParam][index] << std::endl;
      }
      float mean = sum/(float)n;
      float sigma = sqrt( sum2/(float)n - mean*mean);
      borne[iParam].first  = mean - 3*sigma;
      borne[iParam].second = mean + 3*sigma;
      firstRead[iParam] = false;

      initHisto(iParam);

    }
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
};

class FullReader: public ReaderBase{
public:
  FullReader(std::string _filename): ReaderBase(_filename){};

protected:
  virtual void plot(){
    if( nPlotPerCanvas > nVar ) nPlotPerCanvas = nVar;
    for( int iCan = 0; iCan < nVar/nPlotPerCanvas+1*(nVar%nPlotPerCanvas != 0); iCan++){
      TCanvas* can = new TCanvas();
      can -> Divide( nPlotPerCanvas, nPlotPerCanvas );
      for(int iPlot = 0;iPlot < nPlotPerCanvas;iPlot++){
	int iVar = iPlot+iCan*nPlotPerCanvas;
	if( iVar >= nVar ) return;
	for(int jPlot = 0; jPlot <= iPlot; jPlot++){
	  int jVar = jPlot+iCan*nPlotPerCanvas;
	  can -> cd( 1+iPlot*nPlotPerCanvas+jPlot );
                    
	  Stack* st = new Stack(Form("Param #%i",iVar));
	  if( iVar==jVar ){
	    st -> push_back( h[iVar] );
	    //st -> pushVerticalLine( metadata[ Form("realValues_%i",iVar) ] );
	  } else st -> push_back( h2[ std::pair<int, int>(iVar,jVar) ] );
	  st -> draw(gPad);
	}
      }
    }
  }

  virtual void fillHistos(std::string fname, int iParam){
    for(int j = 0;j < metadata["nStep"] ;j++){
      h[iParam] -> Fill( trace[iParam][j] );
      //std::cout << "trace["<<iParam<<"]["<< j << "] = " << trace[iParam][j] << std::endl;
    }

    for(int j = iParam-3;j<iParam && j>=0;j++){
      for(int k = 0;k < metadata["nStep"] ;k++) h2[ std::pair<int,int>(iParam,j) ] -> Fill( trace[iParam][k], trace[j][k] );
    }
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

  virtual void fillHistos(std::string fname, int iParam){
    if(iParam != nVar ) return;

    for(int j = 0;j<nVar;j++){
        for(int k = 0;k < metadata["nStep"] ;k++) h2[ std::pair<int,int>(iParam,j) ] -> Fill( trace[iParam][k], trace[j][k] );
    }

    for(int k = 0;k < metadata["nStep"] ;k++) h[ nVar ] -> Fill( trace[nVar][k] );

  }

  void initHisto(int iParam){
    if(iParam < nVar) return;
    h[iParam] = new TH1F(Form("h_%i",iParam),Form("h_%i",iParam), 100,borne[iParam].first,borne[iParam].second);
        
    for(int j = 0;j<nVar;j++){
      if( firstRead[j] == false ){
        h2[ std::pair<int,int>(iParam,j) ] = new TH2F(Form("h2_%i_%i",iParam,j),Form("h2_%i_%i",iParam,j), 100,borne[iParam].first, borne[iParam].second, 100, borne[j].first, borne[j].second);
      }
    }
  }


};


int main(int argc, char** argv){
  TApplication app("app",&argc,argv,0,-1);

  int c;
  std::string dirname = "test";
  std::cout << "argv[1] : " << argv[1] << std::endl;
  if( argc > 1 ) dirname = argv[1];

  ReaderBase *r = NULL;

  while((c =  getopt(argc, argv, "t:")) != EOF)
    {
      switch (c)
        {
        case 't':
            std::string a(optarg);
            if( a == "full" ){
                r = new FullReader(dirname);
            } else if( a == "lik" ){
                r = new LikelihoodTracer(dirname);
            }
          break;
        }
    }
    
  std::clock_t start = std::clock();

  if( r == NULL ) r = new FullReader(dirname);

  //    r.setMaxSteps(1000);
  r->readAll();

  std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
  app.Run();
  return 0;
}
