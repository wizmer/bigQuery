#ifndef _MCMC__HPP_
#define _MCMC__HPP_

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


#include "generalUtils.hpp"
#include "realisticToyModel.hpp"
#include "ProposalFunction.hpp"



std::ostream& operator<<(std::ostream& os, const SearchSpace& point)
{
    os << "fluxP ";
    for(auto v : point.fluxP) os << v << " ";
    os << "\n";
    os << "fluxD ";
    for(auto v : point.fluxD) os << v << " ";
    return os;
}


template <class Model, class ProposalFunction > 
class MCMC 
{
public:
    MCMC(std::string name = "mcmc.mcmc");
    virtual ~MCMC(){}

    void run(){ loop(); }
    void setSteps(int _nSteps){ nStep = _nSteps; }
    void setVerbose(bool isVerbose){ verbose = isVerbose; }
    void setRegularizationFactor(float _regularizationFactor){
        this -> regularizationFactor = _regularizationFactor;
    }

    float GetLogLikelihood(const SearchSpace & point);


private:

    Model model;

    SearchSpace realValues;
    SearchSpace current_point;
    SearchSpace initialConditions;

    ProposalFunction proposalFunction;

    unsigned int nVar;

    // Raw data 
    std::vector<float> log_likelihood;
    std::vector< std::vector<float> > trace;

    std::string filename;

    float current_log_likelihood;
    float regularizationFactor;

    unsigned int seed;
    unsigned int chunkStepNumber;
    unsigned int nThreads;
    unsigned int nStep;
    unsigned int chunkSize;

    bool verbose;

    static const int maxRAM = 1e9;

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
      
        myfile << "bins " ;
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

    void updateStep(const SearchSpace &proposed_point, const float &proposed_log_likelihood)
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

        current_log_likelihood = GetLogLikelihood(initialConditions);
        
        for( int i = 0; i < nStep; i++)
            {
                if( i%30000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));

                auto proposed_point = proposalFunction.proposePoint(current_point);

                float proposed_log_likelihood = GetLogLikelihood(proposed_point);

                float the_likelihood_ratio = exp(proposed_log_likelihood-current_log_likelihood);

                if( verbose ){

                    std::cout << "fluxP "; for(auto v: proposed_point.fluxP) std::cout << v << " "; std::cout << "\n";
                    std::cout << "fluxD "; for(auto v: proposed_point.fluxD) std::cout << v << " "; std::cout << "\n";

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

#endif
