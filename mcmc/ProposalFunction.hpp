#ifndef _PROPOSALFUNCTION__HPP_
#define _PROPOSALFUNCTION__HPP_
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

std::default_random_engine generator;
std::normal_distribution<float> normalDistrib; // a random gaussian generator

class ProposalFunction {
    std::vector<float> sigmaP, sigmaD;
public: 
    ProposalFunction(const SearchSpace & initialPoint):
        sigmaP(), sigmaD()
    {
        for(auto v : initialPoint.fluxP){
            
            sigmaP.push_back(sqrt(v)>1 ? sqrt(v) : 1);
        }

        for(auto v : initialPoint.fluxD){
            sigmaD.push_back(sqrt(v)>1 ? sqrt(v) : 1 );
        }
    }

    SearchSpace proposePoint(const SearchSpace &previous_point)
    {
        SearchSpace ret = previous_point;
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

#endif
