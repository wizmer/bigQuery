#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <math.h>
#include <cstring>
#include <vector>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "generalUtils.hpp"

std::default_random_engine generator;
std::normal_distribution<float> normalDistrib; // a random gaussian generator

template <class LogP, class ProposalFunction > class MCMC{
public:
    MCMC(std::vector<float> _initialCondition, std::string name = "mcmc.mcmc", std::vector<float> _realValue = std::vector<float>() ){

        
        filename = name;

        if( mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            std::cout << "Cannot create directory: " << filename << std::endl;
            exit(-1);
        }
        
        verbose = false;

        chunkStepNumber = 0;
        nVar = _initialCondition.size();
        nStep=100;
        nThreads = 1;

        chunkSize = maxRAM / ( (nVar+1)*sizeof(float)) / nThreads;

        realValues       = new float[nVar];
        current_point    = new float[nVar];
        initialCondition = new float[nVar];
        log_likelihood   = new float[chunkSize];
        data             = new float[nVar];
        sigma            = new float[nVar];
        trace            = new float*[nVar];

        for(int i = 0;i<nVar;i++){
            trace[i] = new float[chunkSize];
            initialCondition[i] = _initialCondition[i];
            current_point[i] = initialCondition[i];
            realValues[i] = _realValue[i];
            data[i] = _realValue[i];
            sigma[i] = sqrt( _initialCondition[i] );
	}

        // construct a trivial random generator engine from a time-based seed:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(seed);
    }


    void run(){
        loop();
    }

    void setSteps(int _nSteps){
        nStep = _nSteps;
    }

private:
    int chunkSize;
    int nVar;

    std::string filename;

    bool verbose;
    int nStep;

    float current_log_likelihood;

    float *realValues;
    float *sigma;
    float *current_point;
    float *initialCondition;
    float **trace; // the chain containing the log likelihood of all accepted points
    float *log_likelihood; // the chain containing the log likelihood of all accepted points
    float *data;

    // float (*logp)(float*);
    LogP logP;
    ProposalFunction proposalFunction;

    int chunkStepNumber;
    int seed;
    int nThreads;
    static const int maxRAM = 2e9;

    void resetChains(){
        for(int i = 0;i<nVar;i++){
            memset(trace[i],0,chunkSize*sizeof(float));
        }
        memset(log_likelihood,0,chunkSize*sizeof(float)); 
        chunkStepNumber= 0 ;
    }

    void saveMetaData(){
        std::ofstream myfile( filename+"/metadata.txt", std::ios::out);
        myfile << "nStep " << nStep << std::endl;
        myfile << "nVar " << nVar << std::endl;
        myfile << "chunkSize " << chunkSize << std::endl;
        for(int i = 0;i<nVar;i++){
            myfile << "initialPoints_" << i << " " << initialCondition[i] << std::endl;
            myfile << "realValues_" << i << " " << realValues[i] << std::endl;
        }
        myfile.close();
    }

    void saveChunk(){
        static int chunkNumber = 0;
        
        std::cout << "saving chunk" << std::endl;
        for(int i = 0;i<nVar;i++){
            std::stringstream fname;
            fname << filename <<"/par" << i << "_chunk" << chunkNumber << ".bin";
            std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);

            myfile.write((char*)&chunkStepNumber, sizeof(int));
    
            //std::clock_t start = std::clock();

            myfile.write((char*)trace[i], sizeof(float)*chunkStepNumber);
            myfile.close();
	}

        std::stringstream fname;
        fname << filename << "/lik_chunk" << chunkNumber << ".bin";
        std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);
        myfile.write((char*)&chunkStepNumber, sizeof(int));
        myfile.write((char*)log_likelihood, sizeof(float)*chunkStepNumber);
        myfile.close();

        chunkStepNumber = 0;
        chunkNumber++;
    }

    void setSigmas(float *_sigma){
        for(int i = 0;i<nVar;i++){
            sigma[i] = _sigma[i];
	}
    }

    void setVerbose(bool isVerbose){
        verbose = isVerbose;
    }

    void updateStep(float* proposed_point, const float &proposed_log_likelihood){
        if(verbose) std::cout << "accepted" << std::endl;
        for(int i = 0;i<nVar;i++){
            current_point[i] = proposed_point[i];
            current_log_likelihood = proposed_log_likelihood;
        }

        addCurrentPointToChain();
    }

    void addCurrentPointToChain(){
        for(int i = 0;i<nVar;i++){
            trace[i][chunkStepNumber] = current_point[i];
	}
        log_likelihood[chunkStepNumber] = current_log_likelihood;
        chunkStepNumber++;
    }

    void loop(){
        std::cout << "nStep : " << nStep/1e6 << " millions" << std::endl;
        std::cout << "chunkSize : " << chunkSize/1e6 << " MBytes" << std::endl;
        std::cout << "maxRAM : " << maxRAM/1e6 << " MBytes" << std::endl;
        std::cout << "nVar : " << nVar << std::endl;
        std::cout << "sizeof(float) : " << sizeof(float) << std::endl;

        sleep(1);

        saveMetaData();

        current_log_likelihood = LogP::getValue(initialCondition, data, nVar);

        for( int i = 0; i<nStep; i++){
            float proposed_point[nVar];
            if( i%400000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));
            ProposalFunction::proposePoint(current_point, proposed_point, nVar, sigma);
            float proposed_log_likelihood = LogP::getValue(proposed_point, data, nVar);
            float the_likelihood_ratio = exp(proposed_log_likelihood-current_log_likelihood);

            if( verbose ){
                std::cout << "current_log_likelihood : " << current_log_likelihood << std::endl;
                std::cout << "proposed_log_likelihood : " << proposed_log_likelihood << std::endl;
                std::cout << "the_likelihood_ratio : " << the_likelihood_ratio << std::endl;
            }

            // Accept the proposed point if:
            // the likelihood of the proposed point is bigger than previous likelihood
            // OR
            // if smaller, accept it with a chance equal to the likelihood ratio
            // rand() / (float) RAND_MAX, generate a random number between 0 and 1
            if( the_likelihood_ratio > 1 || rand() / (float) RAND_MAX < the_likelihood_ratio ){
                updateStep(proposed_point, proposed_log_likelihood);
            } else addCurrentPointToChain();

            if( chunkStepNumber >= chunkSize) saveChunk();

        }
        
        std::cout << "chunkStepNumber : " << chunkStepNumber << std::endl;
        if( chunkStepNumber > 0 ) saveChunk();
    }
};


 struct DefaultProposalFunction{
     static void proposePoint(float *previous_point, float* proposed_point, const int &nVar, float* sigma){
        for(int i = 0;i<nVar;i++){
            proposed_point[i] = previous_point[i] + sigma[i] * normalDistrib(generator);
        }
    }
};

struct DefaultLogP{
    static float getValue(float *point, float* data, const int &nVar){
        float log = 0;
        for( int i = 0; i<nVar;i++){
            log+= -pow((point[i]-data[i])/sqrt(data[i]),2);
        }
        return log;
    }
};


int main(int argc, char** argv){
    int c;
    int nStep = 0;
    std::string name = "test";

    while((c =  getopt(argc, argv, "n:f:")) != EOF)
        {
            switch (c)
                {
                case 'n':
                    nStep = generalUtils::stringTo<int>(optarg);
                    break;
                case 'f':
                    name = optarg;
                    break;
                }
        }

    std::clock_t start = std::clock();

    int nVar = 62;
    std::vector< float > initialCondition;
    for(int i = 0;i<nVar;i++){
        initialCondition.push_back( 1e6*pow(i+1,-2) );
    }

    MCMC<DefaultLogP, DefaultProposalFunction > a(initialCondition, name, initialCondition);
    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    return 0;
}

