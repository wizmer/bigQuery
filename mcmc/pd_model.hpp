#include <vector>
#include <functional>
#include <iostream>

#include "Matrix.hpp"

class PDModel
{
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    Matrix rgdtF_transposed,  betaF;
    Matrix deltaP, deltaD;
    Matrix observed;
    std::vector<Matrix> matrixBase;
    

    void SetRigidityResolution(const Matrix & matrix);
    void SetBetaResolution    (const Matrix & matrix);

    // Build prediction matrices for all unitary fluxes
    void constuctBaseMatrices();

public:
    static const float mp;
    static const float md;

    // Initializations 
    PDModel( const std::vector<float> & bT, const std::vector<float> & bM, 
             const std::vector<float> & rT, const std::vector<float> & rM,
             const Matrix & betaF, const Matrix & rgdtF);

    static PDModel FromCSVS(const std::string & betaFile, const std::string & rgdtFile, int maxTrueBinNumber = 0 );

    // Getters
    inline std::vector<float> getBetaBinsT(){ return betaBinsT; }
    inline std::vector<float> getBetaBinsM(){ return betaBinsM; }
    inline std::vector<float> getRgdtBinsT(){ return rgdtBinsT; }
    inline std::vector<float> getRgdtBinsM(){ return rgdtBinsM; }
    
    // Predictions
    Matrix GetPrediction( const float* const fluxP,
                          const float* const fluxD  );

    void LoadObservedDataFromFile(const std::string & fname);

    // Perform a linear combination of base matrices
    Matrix GetPredictionFast( const float* const fluxP,
                          const float* const fluxD  );

    float GetLogLikelihood( const float* const fluxP,
                            const float* const fluxD  );

    void GenerateToyObservedData(const std::vector<float> &fluxP,
                                 const std::vector<float> &fluxD  ){
        observed = GetPrediction( &fluxP[0], &fluxD[0]);

        std::cout << "observed.getNrows() : " << observed.getNrows() << std::endl;
        std::cout << "observed.getNcolums() : " << observed.getNcolums() << std::endl;

    }

};

