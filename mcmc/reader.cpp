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

class Reader{
    std::string filename;
public:
    Reader(std::string _filename){
        filename = _filename;
        maxSteps = 0;
        nPlotPerCanvas = 6;
    }

    void setMaxSteps(int _maxSteps){
        maxSteps = _maxSteps;
    }

    void readAll(){
        init();
        // Extract number of variables

        int nChunks = generalUtils::stringTo<int>( generalUtils::exec("ls -lrt " + filename + " | grep -o \"chunk[0-9]*\" | cut -c 6- | sort | uniq | wc -l") );
        std::cout << "nChunks : " << nChunks << std::endl;

        // borne[0] = std::make_pair(999.99e3, 1000.01e3);
        // borne[1] = std::make_pair(200e3, 300e3);
        // borne[2] = std::make_pair(100e3, 100e3);
        // borne[3] = std::make_pair(0, 0);

        for(int jChunk = 0; jChunk<nChunks; jChunk++){
            for(int iVar = 0; iVar<nVar ;iVar++){
                std::string chunkName = Form("%s/par%i_chunk%i.bin", filename.c_str(), iVar, jChunk);
                readBinary( chunkName, iVar );
                for(int j = 0;j < metadata["nStep"] ;j++){
                    h[iVar] -> Fill( trace[iVar][j] );
                    //std::cout << "trace["<<iVar<<"]["<< j << "] = " << trace[iVar][j] << std::endl;
                }

		for(int j = 0;j<3;j++){
		  if(iVar - j< 0) continue;
		  for(int k = 0;k < metadata["nStep"] ;k++) h2[ std::pair<int,int>(iVar,iVar-j) ] -> Fill( trace[iVar][k], trace[iVar-j][k] );
                }
            }
                
        }

        plot();
    }
    
private:
    void plot(){
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
        
    void init(){
        readMetaData();
        trace = new float*[(int)metadata["nVar"]];
        firstRead = new bool[(int)metadata["nVar"]];
        borne = new std::pair<int,int>[nVar];
        for(int i = 0;i<metadata["nVar"];i++){
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
                sum2 += pow(trace[iParam][index], 2);
                //std::cout << "trace["<<iParam<<"]["<< index << "] = " << trace[iParam][index] << std::endl;
            }
            float mean = sum/(float)n;
            float sigma = sqrt( sum2/(float)n - mean*mean);
            borne[iParam].first  = mean - 3*sigma;
            borne[iParam].second = mean + 3*sigma;
            firstRead[iParam] = false;

            h[iParam] = new TH1F(Form("h_%i",iParam),Form("h_%i",iParam), 100,borne[iParam].first,borne[iParam].second);
            
            for(int j = 0;j<nVar;j++){
                if( firstRead[j] == false ){
                    h2[ std::pair<int,int>(iParam,j) ] = new TH2F(Form("h2_%i_%i",iParam,j),Form("h2_%i_%i",iParam,j), 100,borne[iParam].first, borne[iParam].second, 100, borne[j].first, borne[j].second);
                }
            }

        }
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

int main(int argc, char** argv){
    TApplication app("app",&argc,argv,0,-1);

    int c;
    while((c =  getopt(argc, argv, ":a:b:c")) != EOF)
        {
            switch (c)
                {
                case 'a':
                    std::cout << optarg << std::endl;
                    break;
                case 'b':
                    std::cout << optarg << std::endl;
                    break;
                }
        }
    
    std::clock_t start = std::clock();

    std::string dirname = "test";
    if( argc > 1 ) dirname = argv[1];
    Reader r(dirname);
    //    r.setMaxSteps(1000);
    r.readAll();

    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    app.Run();
    return 0;
}
