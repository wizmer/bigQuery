#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <math.h>
#include <cstring>
#include <vector>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "generalUtils.hpp"
#include "pd_model.hpp"

std::default_random_engine generator;
std::normal_distribution<float> normalDistrib; // a random gaussian generator

template <class Model, class ProposalFunction > class MCMC{
public:
    MCMC(std::string name = "mcmc.mcmc" ){
        filename = name;

        if( mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            std::cout << "Cannot create directory: " << filename << std::endl;
            exit(-1);
        }
        
        verbose = true;

        chunkStepNumber = 0;
        initialConditions = model.initialConditions;

        nVar = initialConditions.size();
        nStep=100;
        nThreads = 1;

        chunkSize = maxRAM / ( (nVar+1)*sizeof(float)) / nThreads;


        current_point = initialConditions;
        realValues = model.realValues;

        log_likelihood.reserve(chunkSize);

        for(int i = 0;i<nVar;i++){
            trace.push_back( std::vector<float>(chunkSize) );
            sigma.push_back( sqrt( initialConditions[i] ) );
	}

        // construct a trivial random generator engine from a time-based seed:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(seed);
    }

    void run(){
        loop();
    }

    void setSteps(int _nSteps){
        nStep = _nSteps;
    }

private:
    int chunkSize;
    int nVar;

    std::string filename;

    bool verbose;
    int nStep;

    float current_log_likelihood;

    std::vector<float> realValues;
    std::vector<float> current_point;
    std::vector<float> initialConditions;
    std::vector<float> log_likelihood;
    std::vector<float> data;
    std::vector<float> sigma;
    std::vector< std::vector<float> > trace;

    // float (*logp)(float*);
    Model model;
    ProposalFunction proposalFunction;

    int chunkStepNumber;
    int seed;
    int nThreads;
    static const int maxRAM = 2e9;

    void saveMetaData(){
        std::ofstream myfile( filename+"/metadata.txt", std::ios::out);
        myfile << "nStep " << nStep << std::endl;
        myfile << "nVar " << nVar << std::endl;
        myfile << "chunkSize " << chunkSize << std::endl;
        for(int i = 0;i<nVar;i++){
            myfile << "initialPoints_" << i << " " << initialConditions[i] << std::endl;
            myfile << "realValues_" << i << " " << realValues[i] << std::endl;
        }
        myfile.close();
    }

    void saveChunk(){
        static int chunkNumber = 0;
        
        std::cout << "saving chunk" << std::endl;
        for(int i = 0;i<nVar;i++){
            std::stringstream fname;
            fname << filename <<"/par" << i << "_chunk" << chunkNumber << ".bin";
            std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);

            myfile.write((char*)&chunkStepNumber, sizeof(int));
    
            //std::clock_t start = std::clock();

            myfile.write((char*)&trace[i][0], sizeof(float)*chunkStepNumber);
            myfile.close();
	}

        std::stringstream fname;
        fname << filename << "/lik_chunk" << chunkNumber << ".bin";
        std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);
        myfile.write((char*)&chunkStepNumber, sizeof(int));
        myfile.write((char*)&log_likelihood[0], sizeof(float)*chunkStepNumber);
        myfile.close();

        chunkStepNumber = 0;
        chunkNumber++;
    }

    void setSigmas(float *_sigma){
        for(int i = 0;i<nVar;i++){
            sigma[i] = _sigma[i];
	}
    }

    void setVerbose(bool isVerbose){
        verbose = isVerbose;
    }

    void updateStep(const std::vector<float> &proposed_point, const float &proposed_log_likelihood){
        if(verbose) std::cout << "accepted" << std::endl;
        for(int i = 0;i<nVar;i++){
            current_point[i] = proposed_point[i];
            current_log_likelihood = proposed_log_likelihood;
        }

        addCurrentPointToChain();
    }

    void addCurrentPointToChain(){
        for(int i = 0;i<nVar;i++){
            trace[i][chunkStepNumber] = current_point[i];
	}
        log_likelihood[chunkStepNumber] = current_log_likelihood;
        chunkStepNumber++;
    }

    void loop(){
        std::cout << "nStep : " << nStep/1e6 << " millions" << std::endl;
        std::cout << "chunkSize : " << chunkSize/1e6 << " MBytes" << std::endl;
        std::cout << "maxRAM : " << maxRAM/1e6 << " MBytes" << std::endl;
        std::cout << "nVar : " << nVar << std::endl;
        std::cout << "sizeof(float) : " << sizeof(float) << std::endl;

        sleep(1);

        saveMetaData();

        current_log_likelihood = model.getLogLikelihood(initialConditions, data, nVar);

        for( int i = 0; i<nStep; i++){
            std::vector<float> proposed_point(nVar);
            if( i%400000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));
            ProposalFunction::proposePoint(current_point, proposed_point, nVar, sigma);
            float proposed_log_likelihood = model.getLogLikelihood(proposed_point, data, nVar);
            float the_likelihood_ratio = exp(proposed_log_likelihood-current_log_likelihood);

            if( verbose ){
                std::cout << "current_log_likelihood : " << current_log_likelihood << std::endl;
                std::cout << "proposed_log_likelihood : " << proposed_log_likelihood << std::endl;
                std::cout << "the_likelihood_ratio : " << the_likelihood_ratio << std::endl;
            }

            // Accept the proposed point if:
            // the likelihood of the proposed point is bigger than previous likelihood
            // OR
            // if smaller, accept it with a chance equal to the likelihood ratio
            // rand() / (float) RAND_MAX, generate a random number between 0 and 1
            if( the_likelihood_ratio > 1 || rand() / (float) RAND_MAX < the_likelihood_ratio ){
                updateStep(proposed_point, proposed_log_likelihood);
            } else addCurrentPointToChain();

            if( chunkStepNumber >= chunkSize) saveChunk();

        }
        
        std::cout << "chunkStepNumber : " << chunkStepNumber << std::endl;
        if( chunkStepNumber > 0 ) saveChunk();
    }
};


 struct DefaultProposalFunction{
     static void proposePoint(const std::vector<float> &previous_point, std::vector<float> proposed_point, const int &nVar, const std::vector<float> &sigma){
        for(int i = 0;i<nVar;i++){
            proposed_point[i] = previous_point[i] + sigma[i] * normalDistrib(generator);
        }
    }
};

struct DefaultLogP{
    static float getValue(const std::vector<float> &point, const std::vector<float> &data, const int &nVar){
        float log = 0;
        for( int i = 0; i<nVar;i++){
            log+= -pow((point[i]-data[i])/sqrt(data[i]),2);
        }
        return log;
    }
};


std::vector<float> getBinningFromMatrixFirstRow( std::string matrixFileName ){
    std::string cmd = "cat " + matrixFileName + " | head -n 1 | awk -F',' '{for(i=2;i<=NF;i++) print $i}'";
    return generalUtils::stringTo<float>( generalUtils::splitIntoLines( generalUtils::exec( cmd ) ) );
}

std::vector<float> getBinningFromMatrixFirstColumn( std::string matrixFileName ){
    std::string cmd = "cat " + matrixFileName + " | awk -F',' '{print $1}' | tail -n +2";
    return generalUtils::stringTo<float>( generalUtils::splitIntoLines( generalUtils::exec( cmd ) ) );
}

void fillMatrixFromPandasFile( Matrix &matrix, std::string filename){
    std::vector<std::string> matrixRows = generalUtils::splitIntoLines( generalUtils::exec("cat " + filename + " | tail -n +2 | awk -F',' '{for(i=2;i<NF;i++)printf(\"%s \",$i); printf(\"\\n\")}'") );
    for(int row = 0;row<matrixRows.size()-1;row++){ // Last line is for overflow, hence the -1
        std::cout << "matrixRows[row] : " << matrixRows[row] << std::endl;
        std::vector<float> rowElements = generalUtils::stringTo<float>( generalUtils::split(matrixRows[row], " "));
        for(int column = 0;column<rowElements.size();column++){
            matrix.set(column,row, rowElements[column]);
        }
    }

    matrix.dump();
}


struct RealisticToyModel{

    RealisticToyModel(){
        std::cout << "Hello realistic" << std::endl;
        
        std::string rigidityFile = "datasets/R_resolution.csv";
        std::string betaFile     = "datasets/B_resolution.csv";

        std::vector<float> rgdtBinsT = getBinningFromMatrixFirstRow   ( rigidityFile );
        std::vector<float> rgdtBinsM = getBinningFromMatrixFirstColumn( rigidityFile );

        std::vector<float> betaBinsT = getBinningFromMatrixFirstRow   ( betaFile );
        std::vector<float> betaBinsM = getBinningFromMatrixFirstColumn( betaFile );
    
        model = new PDModel( betaBinsT, betaBinsM, rgdtBinsT, rgdtBinsM);

        Matrix rigidityMatrix( rgdtBinsT.size()-1, rgdtBinsM.size()-1 );
        Matrix betaMatrix    ( betaBinsT.size()-1, betaBinsM.size()-1 );

        fillMatrixFromPandasFile( rigidityMatrix, "datasets/R_resolution.csv");
        fillMatrixFromPandasFile( betaMatrix, "datasets/B_resolution.csv");

        model -> SetRigidityResolution<Matrix>(rigidityMatrix);
        model -> SetBetaResolution<Matrix>(betaMatrix);


        // Set true values of the model
        std::vector<float> toyFluxP, toyFluxD;
        for( int i = 0; i< betaBinsM.size() -1; i++){
            toyFluxP.push_back(100);
            toyFluxD.push_back(100);
        }

        realValues = toyFluxP;
        realValues.insert(realValues.end(), toyFluxD.begin(), toyFluxD.end());

        // Set initial conditions on true value
        initialConditions = realValues;

        // Generate fake data
        model -> GenerateToyObservedData(toyFluxP, toyFluxD);

        
        
    }

    ~RealisticToyModel(){
        delete model;
    }
    

    PDModel *model;
    std::vector<float> initialConditions;
    std::vector<float> realValues;
    
    float getLogLikelihood(const std::vector<float> &point, const std::vector<float> data, const int &nVar){
        float log = model -> GetLogLikelihood( &point[0], &point[initialConditions.size()/2] );

        // # Didn't figure that out yet
        // #firstDerivative = np.diff(np.log(value))
        // #secondDerivative = np.fabs(np.diff(firstDerivative))
        // #smoothness = -(alpha * secondDerivative).sum()
        //         return log #+ smoothness

        return log;
        

    }
};



int main(int argc, char** argv){
    int c;
    int nStep = 0;
    std::string name = "test";

    while((c =  getopt(argc, argv, "n:f:")) != EOF)
        {
            switch (c)
                {
                case 'n':
                    nStep = generalUtils::stringTo<int>(optarg);
                    break;
                case 'f':
                    name = optarg;
                    break;
                }
        }

    std::clock_t start = std::clock();

    MCMC<RealisticToyModel, DefaultProposalFunction > a(name);
    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    return 0;
}

