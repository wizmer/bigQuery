#include "pd_model.hpp"
#include <sstream>
#include <cmath>
#include <chrono>
#include <ctime>

#include "utils/CSV.hpp"
// source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

// Kinematics
const float PDModel::mp = 0.9382;
const float PDModel::md = 1.8756;
float R_from_beta(float beta, float m){ return m*beta/sqrt(1-beta*beta); }
float beta_from_R(float R, float m){ return R/sqrt(R*R+m*m); }


std::vector<float> subBinning(const std::vector<float> &bT, int nBins, int firstBin = 0){
    if( nBins <= 1) return std::vector<float>(bT);

    if( nBins + firstBin > bT.size() ){
        std::cout << "Cannot extract sub binning of size << " << nBins << " starting at: " << firstBin << " using input vector of size: " << bT.size() << std::endl;
        exit(-1);
    }

    std::vector<float> output;
    for(int i = firstBin;i<firstBin+nBins;i++){
        output.push_back( bT[i] );
    }
    return output;
}

//That returns the overlap of the two intervals
float getOverlap(float a0,float a1,float b0,float b1){
    return std::max(0.0F, std::min(a1, b1) - std::max(a0, b0));
}


PDModel::PDModel( 
                 const std::vector<float> & bT, const std::vector<float> & bM,  
                 const std::vector<float> & rT, const std::vector<float> & rM,
                 const Matrix & betaF, const Matrix & rgdtF
                  ): betaBinsT(bT), betaBinsM(bM), betaF(bT.size()-1, bM.size()-1),
                     rgdtBinsT(rT), rgdtBinsM(rM), rgdtF_transposed(rM.size()-1, rT.size()-1),
                     deltaP(bT.size()-1, rT.size()-1), deltaD(bT.size()-1, rT.size()-1),
                     observed(bM.size()-1,rM.size()-1)
{
    for(int bBin=0; bBin < betaBinsT.size() - 1; bBin++)
        {
            for(int rBin=0; rBin < rgdtBinsT.size() - 1; rBin++)
                {
                    float bmin = betaBinsT[bBin], bmax = betaBinsT[bBin+1];
                    float Rmin = rgdtBinsT[rBin], Rmax = rgdtBinsT[rBin+1];
                    float bP1 = beta_from_R(Rmin, mp), bP2 = beta_from_R(Rmax, mp);
                    float bD1 = beta_from_R(Rmin, md), bD2 = beta_from_R(Rmax, md);

                    deltaP.at(bBin, rBin) = getOverlap(bmin,bmax,bP1,bP2)/(bmax-bmin);
                    deltaD.at(bBin, rBin) = getOverlap(bmin,bmax,bD1,bD2)/(bmax-bmin);
                    //if( bBin == betaBinsTs
                }
        }

    SetRigidityResolution(rgdtF);
    SetBetaResolution(betaF);

    constuctBaseMatrices();

}

void PDModel::constuctBaseMatrices(){
    int nVar = betaBinsT.size() -1;

    matrixBase = std::vector<Matrix>(nVar*2);
    for(int i = 0;i<nVar;i++){
        std::vector<float> fakeFluxP(nVar);
        std::vector<float> fakeFluxD(nVar);
        fakeFluxP[i] = 1;
        matrixBase[i] = GetPrediction(&fakeFluxP[0],&fakeFluxD[0]);

        fakeFluxP[i] = 0;
        fakeFluxD[i] = 1;
        matrixBase[i+nVar] = GetPrediction(&fakeFluxP[0],&fakeFluxD[0]);
    }
}

PDModel PDModel::FromCSVS(const std::string & betaFile, const std::string & rgdtFile, int nTrueBins )
{
    std::fstream beta(betaFile), rgdt(rgdtFile);
    std::vector<float> rT, rM, bT, bM;

    Matrix rgdtF = getMatrixAndBins(rgdt, rT, rM).subMatrix(nTrueBins);
    Matrix betaF = getMatrixAndBins(beta, bT, bM).subMatrix(nTrueBins);

    bT = subBinning(bT, nTrueBins+1);
    rT = subBinning(rT, nTrueBins+1);

    PDModel model(bT,bM,rT,rM,betaF,rgdtF);

    return model;
}

void PDModel::SetRigidityResolution(const Matrix & matrix)
{ 
    if( matrix.getNrows() != (rgdtBinsT.size() - 1) ) 
        {
            std::cout <<"Error in " << __FUNCTION__ << "\n";
            std::cout <<"Rigidity resolution matrix size "
                      << "(" << matrix.getNrows()
                      << "," << matrix.getNcolums() << ")\n";
            std::cout <<"is incompatible with rigidity binning size ("
                      << (rgdtBinsT.size() - 1)
                      << ")\n";

            exit(-1);
        }
    rgdtF_transposed = matrix; 
}


void PDModel::SetBetaResolution    (const Matrix & matrix)
{
    if( matrix.getNrows() != (betaBinsT.size() - 1) ) 
        {
            std::cout <<"Error in " << __FUNCTION__ << "\n";
            std::cout <<"Beta resolution matrix size "
                      << "(" << matrix.getNrows()
                      << "," << matrix.getNcolums() << ")\n";
            std::cout <<"is incompatible with beta binning size.\n"
                      << (betaBinsT.size() - 1);
            exit(-1);
        }

    betaF = matrix.Transpose(); 
}


Matrix PDModel::GetPrediction( const float* const fluxP,
                               const float* const fluxD  )
{
    // std::clock_t start = std::clock();
    Matrix fluxMatrixP(deltaP), fluxMatrixD(deltaD);
    
    fluxMatrixP.map([&fluxP](float v, int row, int r){
            return v*fluxP[row];});

    fluxMatrixD.map([&fluxD](float v, int b, int r){return v*fluxD[b];});

    Matrix smearP = betaF.Dot(fluxMatrixP.Dot(rgdtF_transposed));
    Matrix smearD = betaF.Dot(fluxMatrixD.Dot(rgdtF_transposed));

    smearP.map([&smearD](float v, int b, int r){
            //std::cout << v << "\t" << smearD.get(b,r) << std::endl;
            return v + smearD.get(b,r);});

    // std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    return smearP;
}

Matrix PDModel::GetPredictionFast( const float* const fluxP,
                               const float* const fluxD  )
{
    // std::clock_t start = std::clock();
    Matrix output(betaBinsM.size()-1,rgdtBinsM.size()-1);
    int nBinsBetaT = betaBinsT.size()-1;
    for(int i = 0;i<nBinsBetaT;i++){
        float sum = 0;
        output += (matrixBase[i]*fluxP[i]);
        output += (matrixBase[i+nBinsBetaT]*fluxD[i]);
    }
    // std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;

    return output;
}


float PDModel::GetLogLikelihood( const float* const fluxP,
                                 const float* const fluxD  )
{
    Matrix prediction = GetPrediction( fluxP, fluxD );
    // Matrix prediction = GetPredictionFast( fluxP, fluxD );

    // prediction.dump();
    // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << std::endl;
    // observed.dump();
    // exit(-1);

    float ret = prediction.applyAndSum(
                                       [this](float expected , int n, float m){
                                           // std::cout << observed.get(n,m) << "\t" << expected << std::endl;
                                           return observed.get(n,m) * log(expected) - expected;
                                       });

    return ret;
}
