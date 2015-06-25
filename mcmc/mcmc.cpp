#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <math.h>
#include <cstring>

#include <sys/stat.h>
#include <sys/types.h>

class MCMC;

MCMC* gMCMC;

void defaultProposalFunction(float *previous_point, float* proposed_point);
float defaultLogP(float *point);
    
class MCMC{
public:
    MCMC(std::vector<float> _initialCondition, std::string name = "mcmc.mcmc"){
        filename = name;

        if( mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            std::cout << "Cannot create directory: " << filename << std::endl;
            exit(-1);
        }
        
        verbose = false;

        chunkStepNumber = 0;
        nVar = _initialCondition.size();
        nStep=100;

        realValues       = new float[nVar];
        current_point    = new float[nVar];
        initialCondition = new float[nVar];
        sigma            = new float[nVar];
        log_likelihood   = new float[chunkSize];

        trace            = new float*[nVar];

        for(int i = 0;i<nVar;i++){
            trace[i] = new float[chunkSize];
            initialCondition[i] = _initialCondition[i];
            current_point[i] = initialCondition[i];
            sigma[i] = 1;
	}

        nThreads = 1;
        chunkSize = maxRAM / ( (nVar+1)*sizeof(float)) / nThreads;

        std::cout << "chunkSize : " << chunkSize << std::endl;
        std::cout << "maxRAM : " << maxRAM << std::endl;
        std::cout << "nVar : " << nVar << std::endl;
        std::cout << "sizeof(float) : " << sizeof(float) << std::endl;

        // construct a trivial random generator engine from a time-based seed:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(seed);

        gMCMC = this;
        proposalFunction = defaultProposalFunction;
        logp = defaultLogP;
    }

    void defaultProposalMethod(float *previous_point, float* proposed_point){
        for(int i = 0;i<nVar;i++){
            proposed_point[i] = previous_point[i] + sigma[i] * normalDistrib(generator);
        }
    }

    void run(){
        loop();
    }

    float default_logp_method(float *point){
        float log = 0;
        for( int i = 0; i<nVar;i++){
            log+= -pow((point[i]-initialCondition[i])/sqrt(initialCondition[i]),2);
        }
        return log;
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

    float *realValues;
    float *sigma;
    float *current_point;
    float *initialCondition;
    float **trace; // the chain containing the log likelihood of all accepted points
    float *log_likelihood; // the chain containing the log likelihood of all accepted points

    float (*logp)(float*);
    void (*proposalFunction)(float*,float*);

    int chunkStepNumber;
    int seed;
    int nThreads;
    int maxRAM = 2e9;
    std::default_random_engine generator;
    std::normal_distribution<float> normalDistrib; // a random gaussian generator
    

    void resetChains(){
        for(int i = 0;i<nVar;i++){
            memset(trace[i],0,chunkSize*sizeof(float));
        }
        memset(log_likelihood,0,chunkSize*sizeof(float)); 
        chunkStepNumber= 0 ;
    }

    void saveMetaData(){

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

            myfile.write((char*)trace[i], sizeof(float)*chunkStepNumber);
            myfile.close();
	}

        std::stringstream fname;
        fname << filename << "/lik_chunk" << chunkNumber << ".bin";
        std::ofstream myfile( fname.str(), std::ios::out | std::ios::binary);
        myfile.write((char*)&chunkStepNumber, sizeof(int));
        myfile.write((char*)log_likelihood, sizeof(float)*chunkStepNumber);
        myfile.close();

        chunkStepNumber = 0;
        chunkNumber++;
    }

    void setSigmas(float *_sigma){
        for(int i = 0;i<nVar;i++){
            sigma[i] = _sigma[i];
	}
    }

    void setProposalFunction( void(*f)(float*,float*)){
        proposalFunction = f;
    }

    void setLogLikelihoodFunction( float (*f)(float*) ){
        logp = f;
    }


    void setVerbose(bool isVerbose){
        verbose = isVerbose;
    }

    void updateStep(float* proposed_point, const float &proposed_log_likelihood){
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
        // self.metadata={'initialPoints':self.current_point,'nStep':self.nStep,'sigmas':self.sigma,'nVar':self.nVar,'realValues':self.realValues}
        saveMetaData();

        current_log_likelihood = logp(initialCondition);

        for( int i = 0; i<nStep; i++){
            float proposed_point[nVar];
            //if( i%100000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));
            proposalFunction(current_point, proposed_point);
            float proposed_log_likelihood = logp(proposed_point);
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

void defaultProposalFunction(float *previous_point, float* proposed_point){
    gMCMC -> defaultProposalMethod( previous_point, proposed_point );
}

float defaultLogP(float *point){
    gMCMC -> default_logp_method(point);
}


int main(int argc, char** argv){
    std::vector< float > v;
    v.push_back(1000);
    v.push_back(100);
    v.push_back(30);    

    std::string name = "test";
    if( argc > 1 ) name = argv[1];
    MCMC a(v, name);
    a.setSteps(100000000000);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    return 0;
}

