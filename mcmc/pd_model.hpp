#ifndef PD_MODEL_H
#define PD_MODEL_H
#include <vector>
#include <functional>
#include <iostream>

#include "Matrix.hpp"
#include "SearchSpace.hpp"

std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point);

class PDModel
{
    std::vector<MatrixF> matrixBase;
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    MatrixF rgdtF_transposed,  betaF;
    MatrixF deltaP, deltaD;
    MatrixF observed;

    MatrixB mask;

    void SetRigidityResolution(const MatrixF & matrix);
    void SetBetaResolution    (const MatrixF & matrix);
    void SetMask(const MatrixB & _mask);
    void SetMask(const std::string & fname);

    // Build prediction matrices for all unitary fluxes
    void constructBaseMatrices();

    void init(const MatrixF & _betaF, const MatrixF & _rgdtF);
public:
    static const float mp;
    static const float md;

    // Initializations 
    PDModel( const std::vector<float> & bT, const std::vector<float> & bM, 
             const std::vector<float> & rT, const std::vector<float> & rM,
             const MatrixF & betaF, const MatrixF & rgdtF, const MatrixB & mask);

    virtual ~PDModel(){}

    static PDModel FromCSVS(const std::string & betaFile, const std::string & rgdtFile, const std::string & maskFile, int maxTrueBinNumber = 0 );

    // Getters
    inline std::vector<float> getBetaBinsT(){ return betaBinsT; }
    inline std::vector<float> getBetaBinsM(){ return betaBinsM; }
    inline std::vector<float> getRgdtBinsT(){ return rgdtBinsT; }
    inline std::vector<float> getRgdtBinsM(){ return rgdtBinsM; }

    // Predictions
    MatrixF GetPrediction(const SearchSpace & point);

    // Perform a linear combination of base matrices
    MatrixF GetPredictionFast(const SearchSpace & point);

    // Log likelihood
    float GetLogLikelihood(const SearchSpace & point);

    // Regularization term
    float GetRegularizationTerm(const SearchSpace & point);

    // Observed
    void LoadObservedDataFromFile(const std::string & fname);
    void GenerateToyObservedData(const SearchSpace & point){
        observed = GetPrediction(point);
    }
};

std::ostream& operator<<(std::ostream& os, const SearchSpace& point);
#endif //PD_MODEL_H
