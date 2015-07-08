#include "pd_model.hpp"
#include "CSV.hpp"

// source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

// Kinematics
const float PDModel::mp = (float)0.9382;
const float PDModel::md = (float)1.8756;
float R_from_beta(float beta, float m){ return m*beta/(float)sqrt(1-beta*beta); }
float beta_from_R(float R, float m){ return R/(float)sqrt(R*R+m*m); }

std::ostream& operator<<(std::ostream& os, const SearchSpace& point)
{
    os << "fluxP ";
    for(auto v : point.fluxP) os << v << " ";
    os << "\n";
    os << "fluxD ";
    for(auto v : point.fluxD) os << v << " ";
    return os;
}



std::vector< std::pair<float,float> > getLogDerivative(const std::vector< std::pair<float, float> > & point){
    int nVar = point.size();

    std::vector< std::pair<float,float> > derivative(nVar-1);
    for(int i = 0;i<nVar-1;i++){
        derivative[i].first  = ( point[i+1].first - point[i].first )/ (std::log10(point[i+1].second) - std::log10(point[i].second) );
        derivative[i].second = std::sqrt(point[i+1].second * point[i].second);
    }
    return derivative;
}


std::vector<float> subBinning(const std::vector<float> &bT, int nBins, int firstBin = 0){
    if( nBins <= 1) return std::vector<float>(bT);

    if( nBins + firstBin > bT.size() ){
        std::cout << "Cannot extract sub binning of size << " << nBins 
                  << " starting at: " << firstBin 
                  << " using input vector of size: " << bT.size() << std::endl;
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
                 const MatrixF & _betaF, const MatrixF & _rgdtF, const MatrixB & _mask):
    matrixBase((bT.size()-1)*2),
    betaBinsT(bT), betaBinsM(bM), 
    rgdtBinsT(rT), rgdtBinsM(rM),
    rgdtF_transposed(rM.size()-1, rT.size()-1), betaF(bT.size()-1, bM.size()-1),
    deltaP(bT.size()-1, rT.size()-1), deltaD(bT.size()-1, rT.size()-1),
                                                                     observed(bM.size()-1,rM.size()-1), 
                                                                     mask(bM.size()-1,rM.size()-1),
                                                                     regularizationFactor(0),
                                                                     regularizationTerm(0)
{
    SetMask(_mask);
    init(_betaF,_rgdtF);
}

void PDModel::init(const MatrixF & _betaF, const MatrixF & _rgdtF){
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
    
    SetRigidityResolution(_rgdtF);
    SetBetaResolution(_betaF);

    constructBaseMatrices();
}

PDModel PDModel::FromCSVS(const std::string & betaFile, const std::string & rgdtFile, const std::string & maskFile, int nTrueBins )
{
    std::fstream beta(betaFile), rgdt(rgdtFile);
    std::vector<float> rT, rM, bT, bM;

    MatrixF _rgdtF = getMatrixAndBins(rgdt, rT, rM).subMatrix(nTrueBins);
    MatrixF _betaF = getMatrixAndBins(beta, bT, bM).subMatrix(nTrueBins);
    //    MatrixB _mask = getMask(maskFile).subMatrix(nTrueBins);

    MatrixB _mask(bM.size()-1,rM.size()-1,true);
    if( maskFile != "" ) _mask = getMask(maskFile);

    bT = subBinning(bT, nTrueBins+1);
    rT = subBinning(rT, nTrueBins+1);

    PDModel model(bT,bM,rT,rM,_betaF,_rgdtF,_mask);

    return model;
}

void PDModel::SetRigidityResolution(const MatrixF & matrix)
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


void PDModel::SetBetaResolution(const MatrixF & matrix)
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

void PDModel::SetMask(const std::string & maskFile){
    SetMask( getMask(maskFile) );
}

void PDModel::SetMask(const MatrixB & _mask)
{
    if( _mask.getNrows() != (betaBinsM.size() - 1) || _mask.getNcolums() != (rgdtBinsM.size() - 1) ) 
        {
            std::cout <<"Error in " << __FUNCTION__ << "\n";
            std::cout <<"mask matrix size "
                      << "(" << _mask.getNrows()
                      << "," << _mask.getNcolums() << ")\n";
            std::cout <<"is incompatible with measured beta or rgdt binning size"
                      << (betaBinsM.size() - 1) << "," << (rgdtBinsM.size() - 1) << ".\n";
            exit(-1);
        }

    mask = _mask;
}

MatrixF PDModel::GetPrediction(const SearchSpace & point)
{
    // std::clock_t start = std::clock();
    MatrixF fluxMatrixFP(deltaP), fluxMatrixFD(deltaD);
    
    fluxMatrixFP.map([&point](float v, int row){
            return v * point.fluxP[row];});
    fluxMatrixFD.map([&point](float v, int b){
            return v * point.fluxD[b];});

    MatrixF smearP = betaF.Dot(fluxMatrixFP.Dot(rgdtF_transposed));
    MatrixF smearD = betaF.Dot(fluxMatrixFD.Dot(rgdtF_transposed));

    smearP.map([&smearD,this](float v, int b, int r){return (v + smearD.get(b,r))*mask.get(b,r);});

    // std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    return smearP;
}

MatrixF PDModel::GetPredictionFast(const SearchSpace & point)
{
    // std::clock_t start = std::clock();
    MatrixF output(betaBinsM.size()-1,rgdtBinsM.size()-1);
    long unsigned int nBinsBetaT = betaBinsT.size()-1;
    for(int i = 0;i<nBinsBetaT;i++){
        output += (matrixBase[i]*point.fluxP[i]);
        output += (matrixBase[i+nBinsBetaT]*point.fluxD[i]);
    }

    return output;
}

void PDModel::constructBaseMatrices(){
    long unsigned int nVar = rgdtBinsT.size()-1;
    SearchSpace base;
    base.fluxP = std::vector<float>(nVar,0);
    base.fluxD = std::vector<float>(nVar,0);

    for(int i = 0;i<nVar;i++){
        base.fluxP[i] = 1;
        matrixBase[i] = GetPrediction(base);
        base.fluxP[i] = 0;

        base.fluxD[i] = 1;
        matrixBase[i+nVar] = GetPrediction(base);
        base.fluxD[i] = 0;
    }
}


float PDModel::GetLogLikelihood(const SearchSpace & point)
{
    MatrixF prediction = GetPrediction(point);
    float ret = prediction.applyAndSum(
                                       [this](float expected , int n, int m){
                                           if( std::abs(expected) < 1e-99 ) return 0.;
                                           return observed.get(n,m) * log(expected) - expected;
                                       }
                                       );

    ComputeRegularizationTerm(point);
    ret -= regularizationTerm;
    return ret;
}

void PDModel::ComputeRegularizationTerm(const SearchSpace & point){
    regularizationTerm = 0;

    std::vector<std::pair<float,float>> pairFluxEnergyP(point.fluxP.size()), pairFluxEnergyD(point.fluxD.size());
    for(int i = 0;i<pairFluxEnergyP.size();i++)
        pairFluxEnergyP[i] = std::pair<float,float>(point.fluxP[i], std::sqrt(rgdtBinsT[i]*rgdtBinsT[i+1]));

    for(int i = 0;i<pairFluxEnergyD.size();i++)
        pairFluxEnergyD[i] = std::pair<float,float>(point.fluxD[i], std::sqrt(rgdtBinsT[i]*rgdtBinsT[i+1]));


    std::vector< std::pair<float,float> > secondDerivativeP = getLogDerivative( getLogDerivative( pairFluxEnergyP ));
    std::vector< std::pair<float,float> > secondDerivativeD = getLogDerivative( getLogDerivative( pairFluxEnergyD ));

    for(int i = 0;i<secondDerivativeD.size();i++){
        regularizationTerm += std::abs(secondDerivativeP[i].first);
        regularizationTerm += std::abs(secondDerivativeD[i].first);
    }

    regularizationTerm *= regularizationFactor;
}


void PDModel::LoadObservedDataFromFile(const std::string & fname)
{
    std::fstream fs(fname);
    std::vector<std::vector<float> > data;

    for(int b = 0; b < betaBinsM.size(); b++)
        {
            std::vector<float> row;
            for(int r = 0; r < rgdtBinsM.size(); r++)
                {
                    float v = 0;
                    fs >> v;
                    row.push_back(v);
                }
            data.push_back(row);
        }

    MatrixF obs(betaBinsM.size()-1,rgdtBinsM.size()-1);
    obs.Fill(data);
    observed = obs;
}
