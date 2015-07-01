#include <vector>
#include <functional>
#include <iostream>

#include "Matrix.hpp"

class PDModel
{
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    Matrix  rgdtF,  betaF;
    Matrix deltaP, deltaD;
    Matrix observed;
    
public:
    static const float mp;
    static const float md;

    PDModel( const std::vector<float> & bT, const std::vector<float> & bM, 
             const std::vector<float> & rT, const std::vector<float> & rM  );

    Matrix GetPrediction( const float* const fluxP,
                          const float* const fluxD  );

    float GetLogLikelihood( const float* const fluxP,
                            const float* const fluxD  );

    void SetRigidityResolution(const Matrix & matrix);
    void SetBetaResolution    (const Matrix & matrix);
    
    void GenerateToyObservedData(const std::vector<float> &fluxP,
                                 const std::vector<float> &fluxD  ){
        observed = GetPrediction( &fluxP[0], &fluxD[0]);
    }
};

