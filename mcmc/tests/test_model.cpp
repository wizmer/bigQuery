#include "../pd_model.hpp"
#include <cmath>

typedef std::vector<std::vector<float>> V2D;

bool test_model1()
{
    std::vector<float> betaT = { 
        0.5       , 0.62789539, 0.75579077, 0.84989074, 0.91753928,
        0.95511736, 0.9772995 , 0.98817874, 0.99417195, 0.99700241 };

    std::vector<float> rgdtT = {
        0.54167002, 0.75689728, 1.08287816, 1.51314916,  2.16483296,   
        3.02500807, 4.32781998, 6.04743672, 8.65194964, 12.08971682 };
    // 'Magic' bins above without smearing will result in the two diagonals:
    V2D resultData = {
        { 1, 0, 1, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 0, 0, 0, 0, 0 },
        { 0, 0, 1, 0, 1, 0, 0, 0, 0 },
        { 0, 0, 0, 1, 0, 1, 0, 0, 0 },
        { 0, 0, 0, 0, 1, 0, 1, 0, 0 },
        { 0, 0, 0, 0, 0, 1, 0, 1, 0 },
        { 0, 0, 0, 0, 0, 0, 1, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 1, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1 } };
    Matrix R(betaT.size()-1, rgdtT.size()-1);
    R.Fill(resultData);

    // Making a model with two diagonal matrices
    PDModel model(betaT, betaT, rgdtT, rgdtT);

    Matrix unitB(betaT.size()-1, betaT.size()-1);
    Matrix unitR(rgdtT.size()-1, rgdtT.size()-1);

    unitB.map([](float, int n, int m){return n==m?1:0;});
    unitR.map([](float, int n, int m){return n==m?1:0;});

    model.SetRigidityResolution(unitR);
    model.SetBetaResolution    (unitB);

    // Pass unit flux. 
    std::vector<float> unitFlux(betaT.size()-1, 1);
    Matrix prediction = model.GetPrediction(&unitFlux[0],&unitFlux[0]);

    // Test
    
    bool pass = true;
    prediction.map([&pass, &R](float v, int n,int m){
        pass = pass && (std::abs(R.get(n,m) - v) < 0.0001); return v;
    });
    if(!pass)
    {
        std::cout << __FUNCTION__ << " failed:\n";
        prediction.dump(); 
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;
}

bool test_model2()
{
    std::vector<float> betaT = { 
        0.5       , 0.62789539, 0.75579077, 0.84989074, 0.91753928,
        0.95511736, 0.9772995 , 0.98817874, 0.99417195, 0.99700241 };

    std::vector<float> rgdtT = {
        0.54167002, 0.75689728, 1.08287816, 1.51314916,  2.16483296,   
        3.02500807, 4.32781998, 6.04743672, 8.65194964, 12.08971682 };
    // Here I smear down in beta 50/50 
    V2D resultData = {
        { 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0 },
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5 },
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5 },
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5 } };
    Matrix R(betaT.size()-1, rgdtT.size()-1);
    R.Fill(resultData);

    // Making a model with beta smearing 
    PDModel model(betaT, betaT, rgdtT, rgdtT);

    Matrix unitB(betaT.size()-1, betaT.size()-1);
    Matrix unitR(rgdtT.size()-1, rgdtT.size()-1);

    unitB.map([](float, int n, int m){return n==m?0.5:n-m==1?0.5:0;});
    unitR.map([](float, int n, int m){return n==m?1:0;});

    model.SetRigidityResolution(unitR);
    model.SetBetaResolution    (unitB);

    // Pass unit flux. 
    std::vector<float> unitFlux(betaT.size()-1, 1);
    Matrix prediction = model.GetPrediction(&unitFlux[0],&unitFlux[0]);

    // Test
    
    bool pass = true;
    prediction.map([&pass, &R](float v, int n,int m){
        pass = pass && (std::abs(R.get(n,m) - v) < 0.0001); return v;
    });
    if(!pass)
    {
        std::cout << __FUNCTION__ << " failed:\n";
        prediction.dump(); 
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;
}

int main(void)
{
    test_model1();
    test_model2();
    return 0;
}
