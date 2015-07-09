#ifndef _REALISTICTOYMODEL__HPP_
#define _REALISTICTOYMODEL__HPP_

#include <fstream>
#include "generalUtils.hpp"

#include "pd_model.hpp"


// static std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point){
//     long unsigned int nVar = point.size();

//     std::vector< std::pair<float,float> > derivative(nVar-1);
//     for(int i = 0;i<nVar-1;i++){
//         derivative[i].first  = ( point[i+1].first - point[i].first )/ (std::log10(point[i+1].second) - std::log10(point[i].second) );
//         derivative[i].second = std::sqrt(point[i+1].second * point[i].second);
//     }
//     return derivative;
// }

template <typename SearchSpaceType> struct RealisticToyModel: public PDModel<SearchSpaceType>
{
    static const bool isToyMC = true;

    RealisticToyModel():
        PDModel<SearchSpaceType>(PDModel<SearchSpaceType>::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv","datasets/mask.csv" ))
    {
        // Set true values of the model
        for( int i = 0; i < PDModel<SearchSpaceType>::getBetaBinsT().size() - 1; i++){
            this->realValues.fluxP.push_back(10000);
            this->realValues.fluxD.push_back(10000);
        }

        // Generate fake data

        this->GenerateToyObservedData(ModelBase<SearchSpaceType>::realValues);
    }

    virtual void saveMetaData(std::ofstream & myfile){
        myfile << "isToyMC "        << isToyMC << std::endl;
        myfile << "bins " ;
        for(auto v :  this->getBetaBinsT()) myfile << v << " "; 
        myfile << std::endl;
    }
};

template <typename SearchSpaceType> struct RealDataModel: public PDModel<SearchSpaceType>
{
    static const bool isToyMC = false;

    RealDataModel():
        PDModel<SearchSpaceType>(PDModel<SearchSpaceType>::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv", "datasets/mask.csv"))
    {
        // Get initial conditions
        std::string fname = "datasets/initialConditions.txt";
        std::ifstream f(fname);
        if( !f ){
            std::cerr << "No file : " << fname << " found !\nExit !" << std::endl;
            exit(-1);
        }
        while(true){
            float fluxP, fluxD;
            f >> fluxP >> fluxD;
            if( f.eof() ) break;
            this->initialConditions.fluxP.push_back(fluxP);
            this->initialConditions.fluxD.push_back(fluxD);
        }

        if( this->initialConditions.fluxP.size() != this->getBetaBinsT().size()-1 ||
            this->initialConditions.fluxD.size() != this->getBetaBinsT().size()-1 ){
            std::cerr << "Wrong initialConditions in initialConditions.txt" << std::endl;
            std::cerr << "Size do not match: model.getBetaBinsT().size()-1" << std::endl;
            exit(-1);
        }

        // for( int i = 0; i < model.getBetaBinsT().size() - 1; i++){
        //     // putting a realistic power law flux
        //     toyFluxP[i] = 100000* pow(sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]),-2.7);
        //     toyFluxD[i] = toyFluxP[i] / 10.;
        // }

        // Load real data
        this->LoadObservedDataFromFile("datasets/observed_data.txt");
    }

    ~RealDataModel(){}

    virtual void saveMetaData(std::ofstream & myfile){
        myfile << "isToyMC "        << isToyMC << std::endl;
        myfile << "bins " ;
        for(auto v :  this->getBetaBinsT()) myfile << v << " "; 
        myfile << std::endl;
    }
};

template class RealisticToyModel<SearchSpace>;
template class RealDataModel<SearchSpace>;

#endif
