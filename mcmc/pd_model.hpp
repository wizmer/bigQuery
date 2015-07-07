#ifndef PD_MODEL_H
#define PD_MODEL_H
#include <vector>
#include <functional>
#include <iostream>

#include "Matrix.hpp"


class PDModel
{
    std::vector<Matrix> matrixBase;
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    Matrix rgdtF_transposed,  betaF;
    Matrix deltaP, deltaD;
    Matrix observed;
    float regularizationFactor;
    

    void SetRigidityResolution(const Matrix & matrix);
    void SetBetaResolution    (const Matrix & matrix);

    // Build prediction matrices for all unitary fluxes
    //void constuctBaseMatrices();

    // Build prediction matrices for all unitary fluxes
    void constructBaseMatrices();
public:

    struct SearchSpace
    {
        SearchSpace(): fluxP(), fluxD(){}
        ~SearchSpace(){}

        std::vector<float> fluxP;
        std::vector<float> fluxD;
        inline long unsigned int size(){return fluxP.size() + fluxD.size();}
        inline float getRaw(int i){ return i >= fluxP.size()? fluxD[i-fluxP.size()]:fluxP[i]; }

    };

    static const float mp;
    static const float md;

    // Initializations 
    PDModel( const std::vector<float> & bT, const std::vector<float> & bM, 
             const std::vector<float> & rT, const std::vector<float> & rM,
             const Matrix & betaF, const Matrix & rgdtF);

    virtual ~PDModel(){}

    static PDModel FromCSVS(const std::string & betaFile, const std::string & rgdtFile, int maxTrueBinNumber = 0 );

    // Getters
    inline std::vector<float> getBetaBinsT(){ return betaBinsT; }
    inline std::vector<float> getBetaBinsM(){ return betaBinsM; }
    inline std::vector<float> getRgdtBinsT(){ return rgdtBinsT; }
    inline std::vector<float> getRgdtBinsM(){ return rgdtBinsM; }
    
    // Predictions
    Matrix GetPrediction(const SearchSpace & point);

    // Perform a linear combination of base matrices
    Matrix GetPredictionFast(const SearchSpace & point);

    // Log likelihood
    float GetLogLikelihood(const SearchSpace & point);

    // Set regularization factor
    void setRegularization(float _regularizationFactor){
        regularizationFactor = _regularizationFactor;
    }

    // Observed
    void LoadObservedDataFromFile(const std::string & fname);
    void GenerateToyObservedData(const SearchSpace & point){
        observed = GetPrediction(point);
    }

};

std::ostream& operator<<(std::ostream& os, const PDModel::SearchSpace& point);
#endif //PD_MODEL_H
