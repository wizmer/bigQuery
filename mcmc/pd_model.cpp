#include "pd_model.hpp"
#include <sstream>
#include <cmath>

#include "utils/CSV.hpp"
// source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

// Kinematics
const float PDModel::mp = 0.9382;
const float PDModel::md = 1.8756;
float R_from_beta(float beta, float m){ return m*beta/sqrt(1-beta*beta); }
float beta_from_R(float R, float m){ return R/sqrt(R*R+m*m); }

//That returns the overlap of the two intervals
float getOverlap(float a0,float a1,float b0,float b1){
    return std::max(0.0F, std::min(a1, b1) - std::max(a0, b0));
}


PDModel::PDModel( 
    const std::vector<float> & bT, const std::vector<float> & bM,  
    const std::vector<float> & rT, const std::vector<float> & rM
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
        }
    }
}

PDModel PDModel::FromCSVS(const std::string & betaFile, const std::string & rgdtFile )
{
    std::fstream beta(betaFile), rgdt(rgdtFile);
    std::vector<float> rT, rM, bT, bM;

    Matrix rgdtF = getMatrixAndBins(rgdt, rT, rM);
    Matrix betaF = getMatrixAndBins(beta, bT, bM);

    PDModel model(bT,bM,rT,rM);
    model.SetRigidityResolution(rgdtF);
    model.SetBetaResolution(betaF);

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
        std::cout <<"is incompatible with rigidity binning size.\n "
                  << (rgdtBinsT.size() - 1);
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
    Matrix fluxMatrixP(deltaP), fluxMatrixD(deltaD);

    fluxMatrixP.dump();
    std::cout << "\n";
    fluxMatrixD.dump();
    std::cout << "\n";

    fluxMatrixP.map([&fluxP](float v, int b, int r){return v*fluxP[b];});
    fluxMatrixD.map([&fluxD](float v, int b, int r){return v*fluxD[b];});

    Matrix smearP = betaF.Dot(fluxMatrixP.Dot(rgdtF_transposed));
    Matrix smearD = betaF.Dot(fluxMatrixD.Dot(rgdtF_transposed));

    smearP.map([&smearD](float v, int b, int r){return v + smearD.get(b,r);});

    return smearP;
}

float PDModel::GetLogLikelihood( const float* const fluxP,
                                 const float* const fluxD  )
{
    Matrix prediction = GetPrediction( fluxP, fluxD );

    //    prediction.dump();
    // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << std::endl;

    // observed.dump();
    float ret = prediction.applyAndSum(
				  [this](float expected , int n, float m){
				    return observed.get(n,m) * log(expected) - expected;
				  });
    return ret;
}
