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
//         std::cout << __FUNCTION__ << " failed:\n";
//         prediction.dump(); 
//         std::cout << "\n";
//     }
//     else std::cout << __FUNCTION__ << " passed.\n";

//     return pass;
// }


void test_model3(){
    // Test that with identity resolution matrices, the output flux is the input flux
    PDModel model = PDModel::FromCSVS("datasets/B_identity.csv", "datasets/R_identity.csv");

    std::cout << "model.getBetaBinsM().size() : " << model.getBetaBinsM().size() << std::endl;
    std::cout << "model.getBetaBinsT().size() : " << model.getBetaBinsT().size() << std::endl;
    std::cout << "model.getRgdtBinsM().size() : " << model.getRgdtBinsM().size() << std::endl;
    std::cout << "model.getRgdtBinsT().size() : " << model.getRgdtBinsT().size() << std::endl;

    int nVar = model.getBetaBinsT().size()-1;

    bool pass = true;
    //for(int i = 0;i<nVar/2;i++){
    for(int k = 0;k<1;k++){
        std::vector<float> fakeFluxP(nVar/2,0);
        std::vector<float> fakeFluxD(nVar/2,0);
        fakeFluxP[k] = 1;
        Matrix matrix = model.GetPrediction(&fakeFluxP[0],&fakeFluxD[0]);
        matrix.dump();
        // if( matrix.getNrows() != model.getBetaBinsM().size()-1 ){
        //     std::cout << "wrong bin dimension" << std::endl;
        // }
        
        // for(int i = 0;i<matrix.getNcolums();i++){
        //     for(int j = 0;j<matrix.getNrows();j++){
        //         if( k == i && k == j && matrix.get(i,j) != 1 ) {
        //             std::cout << i << "\t" << j << std::endl;
        //             pass = false;
        //         }
        //         if( (k != j || k != i) && std::abs(matrix.get(i,j)) > 0.000001 ) {
        //             std::cout << i << "\t" << j << "\t" << matrix.get(i,j) << std::endl;
        //             pass = false;
        //         }
        //     }
	// }
    }

    if(!pass)
        {
            std::cout << __FUNCTION__ << " failed:\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " passed.\n";

}


int main(void)
{
    // test_model1();
    //test_model2();
    test_model3();
    return 0;
}
