#ifndef PD_MODEL_H
#define PD_MODEL_H
#include <vector>
#include <functional>
#include <iostream>

#include "Matrix.hpp"
//#include "SearchSpace.hpp"
#include "ModelBase.hpp"
#include "SearchSpace.hpp"

std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point);
// std::ostream& operator<<(std::ostream& os, const SearchSpace& point);

template <typename SearchSpaceType> class PDModel: public ModelBase<SearchSpaceType>
{
    std::vector<MatrixF> matrixBase;
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    MatrixF rgdtF_transposed,  betaF;
    MatrixF deltaP, deltaD;
    MatrixF observed;

    MatrixB mask;

    float regularizationFactor;
    float regularizationTerm;


    void SetRigidityResolution(const MatrixF & matrix);
    void SetBetaResolution    (const MatrixF & matrix);

    // Build prediction matrices for all unitary fluxes
    void constructBaseMatrices();

    void init(const MatrixF & _betaF, const MatrixF & _rgdtF);
public:
    static constexpr float mp = (float)0.9382;
    static constexpr float md = (float)1.8756;

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
    float regularizationTermAccessor(){ return regularizationTerm; }

    // Predictions
    MatrixF GetPrediction(const SearchSpaceType & point);
    
    // Perform a linear combination of base matrices
    MatrixF GetPredictionFast(const SearchSpaceType & point);
    
    // Log likelihood
    virtual float GetLogLikelihood(const SearchSpaceType & point);
    
    // Regularization term
    void ComputeRegularizationTerm(const SearchSpaceType & point);

    
    void setRegularizationFactor(float _regularizationFactor){
        this -> regularizationFactor = _regularizationFactor;
    }
    
    // Observed
    void LoadObservedDataFromFile(const std::string & fname);
    void GenerateToyObservedData(const SearchSpaceType & point){
        observed = GetPrediction(point);
    }

    void SetMask(const std::string & maskFile);
    void SetMask(const MatrixB & _mask);

};


#endif //PD_MODEL_H

