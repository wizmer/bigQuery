#ifndef _REALISTICTOYMODEL__HPP_
#define _REALISTICTOYMODEL__HPP_

#include <fstream>
#include "generalUtils.hpp"

#include "pd_model.hpp"

static std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point){
    int nVar = point.size();

    std::vector< std::pair<float,float> > derivative(nVar-1);
    for(int i = 0;i<nVar-1;i++){
        derivative[i].first  = ( point[i+1].first - point[i].first )/ (std::log10(point[i+1].second) - std::log10(point[i].second) );
        derivative[i].second = std::sqrt(point[i+1].second * point[i].second);
    }
    return derivative;
}

struct RealisticToyModel: public PDModel
{
    SearchSpace initialConditions;
    SearchSpace realValues;
    static const bool isToyMC = true;

    RealisticToyModel():
        PDModel(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv"))
    {
        // Set true values of the model
        for( int i = 0; i < getBetaBinsT().size() - 1; i++){
            realValues.fluxP.push_back(10000);
            realValues.fluxD.push_back(10000);
        }


        // Set initial conditions on true value
        initialConditions = realValues;

        // Generate fake data

        GenerateToyObservedData(realValues);
    }
};

struct RealDataModel: public PDModel
{
    SearchSpace initialConditions;
    SearchSpace realValues;
    static const bool isToyMC = false;

    RealDataModel():
        PDModel(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv"))
    {
        // Get initial conditions
        std::ifstream f("datasets/initialConditions.txt");
        while(true){
            float fluxP, fluxD;
            f >> fluxP >> fluxD;
            if( f.eof() ) break;
            initialConditions.fluxP.push_back(fluxP);
            initialConditions.fluxD.push_back(fluxD);
        }

        if( initialConditions.fluxP.size() != getBetaBinsT().size()-1 || initialConditions.fluxD.size() != getBetaBinsT().size()-1 ){
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
        LoadObservedDataFromFile("datasets/observed_data.txt");
    }
};

#endif
