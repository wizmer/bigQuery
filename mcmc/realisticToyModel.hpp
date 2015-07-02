#include "pd_model.hpp"

struct RealisticToyModel: public PDModel
{
    SearchSpace initialConditions;
    SearchSpace realValues;
    static const bool isToyMC = true;

    RealisticToyModel():
        PDModel(PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv",16))
    {
        // Set true values of the model
        for( int i = 0; i < model.getBetaBinsT().size() - 1; i++){
            realValues.fluxP.push_back(10000);
            realValues.fluxD.push_back(10000);
        }

        // Set initial conditions on true value
        initialConditions = realValues;

        // Generate fake data
        model.GenerateToyObservedData(realValues);
    }
};

