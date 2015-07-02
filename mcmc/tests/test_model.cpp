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


    Matrix unitB(betaT.size()-1, betaT.size()-1);
    Matrix unitR(rgdtT.size()-1, rgdtT.size()-1);

    unitB.map([](float, int n, int m){return n==m?1:0;});
    unitR.map([](float, int n, int m){return n==m?1:0;});

    // Making a model with two diagonal matrices
    PDModel model(betaT, betaT, rgdtT, rgdtT, unitB, unitR);

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
        std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m:\n";
        prediction.dump(); 
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

    return pass;
}

// bool test_model2()
// {
//     std::vector<float> betaT = { 
//         0.5       , 0.62789539, 0.75579077, 0.84989074, 0.91753928,
//         0.95511736, 0.9772995 , 0.98817874, 0.99417195, 0.99700241 };

//     std::vector<float> rgdtT = {
//         0.54167002, 0.75689728, 1.08287816, 1.51314916,  2.16483296,   
//         3.02500807, 4.32781998, 6.04743672, 8.65194964, 12.08971682 };
//     // Here I smear down in beta 50/50 
//     V2D resultData = {
//         { 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//         { 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 },
//         { 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0 },
//         { 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0 },
//         { 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0 },
//         { 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0 },
//         { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5 },
//         { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5 },
//         { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5 } };
//     Matrix R(betaT.size()-1, rgdtT.size()-1);
//     R.Fill(resultData);

//     // Making a model with beta smearing 
//     PDModel model(betaT, betaT, rgdtT, rgdtT);

//     Matrix unitB(betaT.size()-1, betaT.size()-1);
//     Matrix unitR(rgdtT.size()-1, rgdtT.size()-1);

//     unitB.map([](float, int n, int m){return n==m?0.5:n-m==1?0.5:0;});
//     unitR.map([](float, int n, int m){return n==m?1:0;});

//     model.SetRigidityResolution(unitR);
//     model.SetBetaResolution    (unitB);

//     // Pass unit flux. 
//     std::vector<float> unitFlux(betaT.size()-1, 1);
//     Matrix prediction = model.GetPrediction(&unitFlux[0],&unitFlux[0]);

//     // Test
    
//     bool pass = true;
//     prediction.map([&pass, &R](float v, int n,int m){
//         pass = pass && (std::abs(R.get(n,m) - v) < 0.0001); return v;
//     });
//     if(!pass)
//     {
//         std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m:\n";
//         prediction.dump(); 
//         std::cout << "\n";
//     }
//     else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

//     return pass;
// }


void test_model3(){
    // Test that with identity resolution matrices, the output flux is the input flux
    PDModel model = PDModel::FromCSVS("datasets/B_identity.csv", "datasets/R_identity.csv",10);

    int nVar = model.getBetaBinsT().size()-1;

    bool pass = true;
    //for(int k = 0;k<nVar;k++){
    for(int k = 0;k<1;k++){
        std::vector<float> fakeFluxP(nVar,0);
        std::vector<float> fakeFluxD(nVar,0);
 
        fakeFluxP[k] = 1;
        Matrix matrix = model.GetPrediction(&fakeFluxP[0],&fakeFluxD[0]);
  
        if( matrix.getNrows() != model.getBetaBinsM().size()-1 ){
            std::cout << "wrong beta bin dimension" << std::endl;
        }

        if( matrix.getNcolums() != model.getRgdtBinsM().size()-1 ){
            std::cout << "wrong rigidity bin dimension" << std::endl;
        }
        
        for(int i = 0;i<matrix.getNrows();i++){
            for(int j = 0;j<matrix.getNcolums();j++){
                if( k == i && k == j && std::abs(matrix.get(i,j)-1) > 0.002 ) {
                    std::cout << i << "\t" << j << "\t" << matrix.get(i,j) << std::endl;
                    pass = false;
                }
                if( (k != j || k != i) && std::abs(matrix.get(i,j)) > 0.002 ) {
                    std::cout << i << "\t" << j << "\t" << matrix.get(i,j) << std::endl;
                    pass = false;
                }
            }
	}
    }

    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m:\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

}


// Test that GetPredictionFast works correctly
void test_model4(){
    PDModel model = PDModel::FromCSVS("datasets/B_resolution.csv", "datasets/R_resolution.csv",10);
    int nVar = model.getBetaBinsT().size()-1;

    bool pass = true;
    std::vector<float> fakeFluxP(nVar,0);
    std::vector<float> fakeFluxD(nVar,0);
    for(int i = 0;i<nVar;i++){
        fakeFluxP[i] = i;
        fakeFluxD[i] = i;
    }

    Matrix A = model.GetPrediction(&fakeFluxP[0],&fakeFluxD[0]);
    Matrix B = model.GetPredictionFast(&fakeFluxP[0],&fakeFluxD[0]);

    for(int i = 0;i<A.getNrows();i++){
        for(int j = 0;j<A.getNcolums();j++){
            if( std::abs(A.get(i,j) - B.get(i,j)) > 0.000001 ){
                std::cout << i << "\t" << j << "\t" << A.get(i,j) << "\t" << B.get(i,j) << std::endl;
                pass = false;
            }
        }
    }

    
    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m:\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";
    
}

int main(void)
{
    test_model1();
    //test_model2();
    test_model3();
    test_model4();
    return 0;
}
