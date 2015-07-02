#ifndef _REALISTICTOYMODEL__HPP_
#define _REALISTICTOYMODEL__HPP_

#include "generalUtils.hpp"

struct RealisticToyModel
{
    PDModel model;
    std::vector<float> initialConditions;
    std::vector<float> realValues;
    static const bool isToyMC = true;
    float alpha;

    RealisticToyModel():
        model(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv",13)), alpha(0)
    {
        // Set true values of the model
        std::vector<float> toyFluxP, toyFluxD;
        for( int i = 0; i < model.getBetaBinsT().size() - 1; i++){
            // putting a realistic power law flux
            toyFluxP.push_back( 100000* pow(sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]),-2.7));
            toyFluxD.push_back( toyFluxP[i] / 10. );
        }

        realValues = toyFluxP;
        realValues.insert(realValues.end(), toyFluxD.begin(), toyFluxD.end());

        // Set initial conditions on true value
        initialConditions = realValues;

        // Generate fake data
        model.GenerateToyObservedData(toyFluxP, toyFluxD);
    }

    float getLogLikelihood(const std::vector<float> &point, const std::vector<float> data, const int &nVar){
        // for(int i = 0;i<point.size();i++){
        //std::cout << "point["<<i<<"] : " << point[i] << std::endl;
        // }
        float log = model.GetLogLikelihood( &point[0], &point[initialConditions.size()/2] );

        std::vector<std::pair<float,float>> fluxP(point.size()/2), fluxD(point.size()/2);
        for(int i = 0;i<point.size()/2;i++){
            fluxP[i] = std::pair<float,float>(point[i]               , std::sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]));
            fluxD[i] = std::pair<float,float>(point[i+point.size()/2], std::sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]));
	}

        std::vector< std::pair<float,float> > secondDerivativeP = getLogDerivative( getLogDerivative( fluxP ));
        std::vector< std::pair<float,float> > secondDerivativeD = getLogDerivative( getLogDerivative( fluxD ));
        float smoothness = 0;
        for(int i = 0;i<secondDerivativeD.size();i++){
            smoothness += std::abs(secondDerivativeP[i].first);
            smoothness += std::abs(secondDerivativeD[i].first);
	}
        //std::cout << log << "\t" << smoothness << std::endl;
        return log - alpha * smoothness;
    }

    static std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point){
        int nVar = point.size();

        std::vector< std::pair<float,float> > derivative(nVar-1);
        for(int i = 0;i<nVar-1;i++){
            derivative[i].first  = ( point[i+1].first - point[i].first )/ (std::log10(point[i+1].second) - std::log10(point[i].second) );
            derivative[i].second = std::sqrt(point[i+1].second * point[i].second);
        }
        return derivative;
    }

    void setRegularization(float _alpha){
        alpha = _alpha;
    }

};

struct RealDataModel
{
    PDModel model;
    std::vector<float> initialConditions;
    std::vector<float> realValues;
    static const bool isToyMC = false;
    float alpha;

    RealDataModel():
        model(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv",13)), alpha(0)
    {
        // Set true values of the model
        std::vector<float> toyFluxP(model.getBetaBinsT().size() - 1);
        toyFluxP = {
            56364.9,
            48061.1,
            40230.2,
            32367.6,
            24068.8,
            16239.9,
            10361  ,
            7038.97,
            5581.88,
            4386.49,
            3266.27,
            2134.53,
            949.041};

        std::vector<float> toyFluxD(model.getBetaBinsT().size() - 1);
        toyFluxD = {
            6820.49,
            5891.42,
            4951.75,
            3971.72,
            3043.72,
            2330.2 ,
            1824.86,
            1321.47,
            815.991,
            424.332,
            195.687,
            92.8731,
            22.2737};

    // for( int i = 0; i < model.getBetaBinsT().size() - 1; i++){
    //         // putting a realistic power law flux
    //         toyFluxP.push_back( 100000* pow(sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]),-2.7));
    //         toyFluxD.push_back( toyFluxP[i] / 10. );
    // }

        initialConditions = toyFluxP;
        initialConditions.insert(initialConditions.end(), toyFluxD.begin(), toyFluxD.end());

        // Load real data
        model.LoadObservedDataFromFile("datasets/observed_data.txt");
    }

    float getLogLikelihood(const std::vector<float> &point, const std::vector<float> data, const int &nVar){
        // for(int i = 0;i<point.size();i++){
        //std::cout << "point["<<i<<"] : " << point[i] << std::endl;
        // }
        float log = model.GetLogLikelihood( &point[0], &point[initialConditions.size()/2] );

        std::vector<std::pair<float,float>> fluxP(point.size()/2), fluxD(point.size()/2);
        for(int i = 0;i<point.size()/2;i++){
            fluxP[i] = std::pair<float,float>(point[i]               , std::sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]));
            fluxD[i] = std::pair<float,float>(point[i+point.size()/2], std::sqrt(model.getRgdtBinsT()[i]*model.getRgdtBinsT()[i+1]));
	}

        std::vector< std::pair<float,float> > secondDerivativeP = getLogDerivative( getLogDerivative( fluxP ));
        std::vector< std::pair<float,float> > secondDerivativeD = getLogDerivative( getLogDerivative( fluxD ));
        float smoothness = 0;
        for(int i = 0;i<secondDerivativeD.size();i++){
            smoothness += std::abs(secondDerivativeP[i].first);
            smoothness += std::abs(secondDerivativeD[i].first);
	}
        //std::cout << log << "\t" << smoothness << std::endl;
        return log - alpha * smoothness;
    }

    static std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point){
        int nVar = point.size();

        std::vector< std::pair<float,float> > derivative(nVar-1);
        for(int i = 0;i<nVar-1;i++){
            derivative[i].first  = ( point[i+1].first - point[i].first )/ (std::log10(point[i+1].second) - std::log10(point[i].second) );
            derivative[i].second = std::sqrt(point[i+1].second * point[i].second);
        }
        return derivative;
    }

    void setRegularization(float _alpha){
        alpha = _alpha;
    }

};

#endif
