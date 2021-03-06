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

struct RealisticToyModel: public PDModel
{
    static const bool isToyMC = true;

    RealisticToyModel():
       PDModel(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv","datasets/mask.csv" ))
    {
        // Set true values of the model
        for( int i = 0; i < getBetaBinsT().size() - 1; i++){
            realValues.fluxP.push_back(10000);
            realValues.fluxD.push_back(10000);
        }

        // Generate fake data

        GenerateToyObservedData(realValues);
    }

    virtual void saveMetaData(std::ofstream & myfile){
        myfile << "isToyMC "        << isToyMC << std::endl;
        myfile << "bins " ;
        for(auto v :  getBetaBinsT()) myfile << v << " "; 
        myfile << std::endl;
    }

};

struct RealDataModel: public PDModel
{
    static const bool isToyMC = false;

    RealDataModel():
        PDModel(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv", "datasets/mask.csv"))
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

    ~RealDataModel(){}

    virtual void saveMetaData(std::ofstream & myfile){
        myfile << "isToyMC "        << isToyMC << std::endl;
        myfile << "bins " ;
        for(auto v :  getBetaBinsT()) myfile << v << " "; 
        myfile << std::endl;
    }
};

#endif
