#include <ctime>
#include <iostream>
#include <iosfwd>
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
#include "realisticToyModel.hpp"

std::default_random_engine generator;
std::normal_distribution<float> normalDistrib; // a random gaussian generator

class ProposalFunction {
    std::vector<float> sigmaP, sigmaD;
public: 
    typedef PD_Model::SearchSpace Space;
    ProposalFunction(const Space & initialPoint)
    {
        for(auto v : initialPoint.fluxP) sigmaP.push_back(v);
        for(auto v : initialPoint.fluxD) sigmaD.push_back(v);
    }

    Space proposePoint(const Space &previous_point)
    {
        Space ret = previous_point;
        for(int i = 0; i < sigmaP.size(); i++) {
            do { ret.fluxP[i] = previous_point.fluxP[i] + sigmaP[i] * normalDistrib(generator); } 
            while( ret.fluxP[i] < 0);
        }
        for(int i = 0; i < sigmaD.size(); i++) {
            do { ret.fluxD[i] = previous_point.fluxD[i] + sigmaD[i] * normalDistrib(generator); } 
            while( ret.fluxD[i] < 0);
        }
        return ret;
    }
};

ostream& operator<<(ostream& os, const PD_Model::SearchSpace& point)
{
    os << "fluxP ";
    for(v : point.fluxP) os << v << " ";
    os << "\n";
    os << "fluxD ";
    for(v : point.fluxD) os << v << " ";
}



template <class Model, class ProposalFunction > 
class MCMC 
{
public:
    MCMC(std::string name = "mcmc.mcmc"):
        initialConditions(model.initialConditions),
        realValues(model.realValues),
        proposalFunction(model.initialConditions)
    {
        filename = name;

        if( mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            int i = 0;
            while( mkdir( (filename+generalUtils::toString(i)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0 ){
                i++;
            }
            filename = filename+generalUtils::toString(i);
        }
        
        verbose = false;

        chunkStepNumber = 0;
        nVar = initialConditions.size();
        nStep=100;
        nThreads = 1;
        chunkSize = maxRAM / ( (nVar+1)*sizeof(float)) / nThreads;

        current_point = initialConditions;
        log_likelihood.reserve(chunkSize);

        for(int i = 0; i < nVar; i++)
            trace.push_back( std::vector<float>(chunkSize) );

        // construct a trivial random generator engine from a time-based seed:
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(seed);

        // Save the predicted (BetaMeas,RMeas) matrix for true flux

        //(model.model -> GetPrediction(&model.realValues[0],&model.realValues[nVar/2])).save(filename+"/predictedMatrix");
    }

    void run(){ loop(); }
    void setSteps(int _nSteps){ nStep = _nSteps; }
    void setVerbose(bool isVerbose){ verbose = isVerbose; }

private:
    int chunkSize;
    int nVar;

    std::string filename;

    bool verbose;
    int nStep;
    unsigned seed;
    float current_log_likelihood;

    Model model;

    Model::SearchSpace realValues;
    Model::SearchSpace current_point;
    Model::SearchSpace initialConditions;

    ProposalFunction proposalFunction;

    // Raw data 
    std::vector<float> log_likelihood;
    std::vector< std::vector<float> > trace;

    int chunkStepNumber;
    int nThreads;

    static const int maxRAM = 2e9;

    void saveMetaData()
    {
        std::ofstream myfile( filename+"/metadata.txt", std::ios::out);
        myfile << "nStep "          << nStep         << std::endl;
        myfile << "nVar "           << nVar          << std::endl;
        myfile << "chunkSize "      << chunkSize     << std::endl;
        myfile << "seed "           << seed          << std::endl;
        myfile << "isToyMC "        << model.isToyMC << std::endl;
        myfile << "initialPoint"    << std::endl;
        myfile << initialConditions << std::endl;
        myfile << "realValues"      << std::endl;
        myfile << realValues        << std::endl;
      
        myfile << "bins " 
        for(auto v :  model.getBetaBinsT()) myfile << v << " "; 
        myfile << std::endl;


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

            myfile.write((char*)&trace[i][0], sizeof(float)*chunkStepNumber);
            myfile.close();
        }

        std::stringstream fname;
        fname << filename << "/lik_chunk" << chunkNumber << ".bin";
        std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);
        myfile.write((char*)&chunkStepNumber, sizeof(int));
        myfile.write((char*)&log_likelihood[0], sizeof(float)*chunkStepNumber);
        myfile.close();

        chunkStepNumber = 0;
        chunkNumber++;
    }

    void updateStep(const Model::SearchSpace &proposed_point, const float &proposed_log_likelihood)
    {
        if(verbose) std::cout << "accepted" << std::endl;

        current_point          = proposed_point;
        current_log_likelihood = proposed_log_likelihood;

        addCurrentPointToChain();
    }

    void addCurrentPointToChain()
    {
        log_likelihood[chunkStepNumber] = current_log_likelihood;
        for(int i = 0; i < nVar; i++)
            trace[i][chunkStepNumber] = current_point.getRaw(i);
        
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

        current_log_likelihood = model.GetLogLikelihood(initialConditions);

        for( int i = 0; i < nStep; i++)
        {
            if( i%10000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));

            auto proposed_point = proposalFunction.proposePoint(current_point);

            float proposed_log_likelihood = model.GetLogLikelihood(proposed_point);
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


int main(int argc, char** argv){
    int c;
    int nStep = 0;
    std::string name = "test";
    bool verbose = false;

    while((c =  getopt(argc, argv, "n:f:v")) != EOF)
        {
            switch (c)
                {
                case 'n':
                    nStep = generalUtils::stringTo<int>(optarg);
                    break;
                case 'f':
                    name = optarg;
                    break;
                case 'v':
                    verbose = true;
                    break;
                }
        }

    std::clock_t start = std::clock();

    MCMC<RealisticToyModel, ProposalFunction > a(name);
    a.setVerbose(verbose);
    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    std::cout << "Processing of : " << name << " done !" << std::endl;
    return 0;
}

