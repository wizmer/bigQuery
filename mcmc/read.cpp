#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "TH1F.h"
#include "TApplication.h"

TH1F* readBinary(std::string fname){
    std::cout << "binary" << std::endl;
    std::ifstream f(fname.c_str(), std::ios::binary );

    float sum = 0;

    unsigned int N;

    f.read((char*)&N, sizeof(unsigned int));

    std::vector<float> a(N);
    std::clock_t start = std::clock();
    f.read((char*)&a[0], sizeof(float)*N);

    TH1F* h = new TH1F("h","h",100,0,0);

    for(int i = 0;i<N;i++){
        std::cout << "a[i] : " << a[i] << std::endl;
        h -> Fill( a[i] );
    }

    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    std::cout << "sum : " << sum << std::endl;
    return h;
}

int main(int argc, char** argv){
    TApplication app("app",&argc,argv);
    TH1F* h = readBinary("exemple_0.bin");
    h -> Draw();
    app.Run();
    return 0;
}
