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
        readMetaData();

        trace = new float*[nVar];
        for(int i = 0;i<metadata["nVar"];i++){
            trace[i] = new float[chunkSize];
	}


        readAll();
    }

private:
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
        
    }
    
    std::vector<float> readBinary(std::string fname, int iParam){
        std::cout << "binary" << std::endl;
        std::ifstream f(fname.c_str(), std::ios::binary );

        float sum = 0;

        unsigned int N;

        f.read((char*)&N, sizeof(unsigned int));
        f.read((char*)trace[iParam], sizeof(float)*N);
    }

    void readAll(){
        // Extract number of variables

        int nChunks = generalUtils::stringTo<int>( generalUtils::exec("ls -lrt " + filename + " | grep -o \"chunk[0-9]*\" | cut -c 6- | sort | uniq | wc -l") );
        std::cout << "nChunks : " << nChunks << std::endl;

        std::map< int, TH1* > h;
        for(int i = 0;i<metadata["nVar"];i++){
            h[i] = new TH1F(Form("h_%i",i),Form("h_%i",i), 100,0,0);
        }

        delete h[0];
        h[0] = new TH1F(Form("h_%i",0),Form("h_%i",0), 100,800,1200);

        for(int jChunk = 0; jChunk<nChunks; jChunk++){
            for(int iVar = 0;iVar<metadata["nVar"];iVar++){
                    std::string chunkName = Form("%s/par%i_chunk%i.bin", filename.c_str(), iVar, jChunk);
                    readBinary( chunkName, h[iVar] );
                }
            }
        }

        TCanvas* can = new TCanvas();
        can -> Divide( metadata["nVar"], metadata["nVar"] );
        for(int i = 0;i<metadata["nVar"];i++){
            for(int j = 0; j < metadata["nVar"]; j++){
                can -> cd( 1+i+metadata["nVar"]+j );
                Stack* st = new Stack(Form("Param #%i",i));
                if( i==j ) st -> push_back( h[i] );
                else st -> push_back( h2[i] );
                st -> pushVerticalLine( metadata[ Form("realValues_%i",i) ] );
                st -> draw();
            }
        }
    }
};

int main(int argc, char** argv){
    TApplication app("app",&argc,argv,0,-1);
    std::clock_t start = std::clock();

    std::string dirname = "test";
    if( argc > 1 ) dirname = argv[1];

    Reader r(dirname);
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    app.Run();
    return 0;
}
