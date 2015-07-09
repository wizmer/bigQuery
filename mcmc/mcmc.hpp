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

#include "utils/generalUtils.hpp"
#include "utils/SearchSpace.hpp"
#include "realisticToyModel.hpp"
#include "ProposalFunction.hpp"


template <class ProposalFunction, class SearchSpaceType> 
class MCMC 
{
public:
    MCMC(std::string name, ModelBase<SearchSpaceType> & _model);
    virtual ~MCMC(){}

    void run(){ loop(); }
    void setSteps(int _nSteps){ nStep = _nSteps; }
    void setVerbose(bool isVerbose){ verbose = isVerbose; }
    //    void setSpectator(std::function<float(const ModelBase&)>);

private:

ModelBase<SearchSpaceType> & model;
    
SearchSpaceType realValues;
    SearchSpaceType current_point;
    SearchSpaceType initialConditions;

    ProposalFunction proposalFunction;

    unsigned int nVar;
    unsigned int nThreads;
unsigned int chunkSize;

// Raw data 
// std::vector<std::function<float(const ModelBase<SearchSpaceBase>&)>> spectators;
std::vector<float> log_likelihood;
    std::vector< std::vector<float> > trace;

    std::string filename;

    float current_log_likelihood;

    unsigned int seed;
    unsigned int chunkStepNumber;
    unsigned int nStep;

    bool verbose;

    static const int maxRAM = 1e9;

    void saveMetaData()
    {
        std::ofstream myfile( filename+"/metadata.txt", std::ios::out);

        myfile << "nStep "          << nStep         << std::endl;
        myfile << "nVar "           << nVar          << std::endl;
        myfile << "chunkSize "      << chunkSize     << std::endl;
        myfile << "seed "           << seed          << std::endl;
        // myfile << "initialPoint"    << std::endl;
        // myfile << initialConditions << std::endl;
        // myfile << "realValues"      << std::endl;
        // myfile << realValues        << std::endl;
      
        model.saveMetaData(myfile);

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

    void updateStep(const SearchSpaceType &proposed_point, const float &proposed_log_likelihood)
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

    void loop();
};

#endif

