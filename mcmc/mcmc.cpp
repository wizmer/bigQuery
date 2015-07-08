#include "mcmc.hpp"

template <class Model, class ProposalFunction > MCMC<Model, ProposalFunction >::MCMC(std::string name, Model _model):
    model(_model),
    realValues(model.realValues),
    current_point(model.initialConditions),
    initialConditions(model.initialConditions),
    proposalFunction(model.initialConditions),
    nVar(initialConditions.size()),
    nThreads(1),
    chunkSize(maxRAM / ( (nVar+1)*sizeof(float)) / nThreads),
    log_likelihood(chunkSize),
    trace(nVar),
    filename(name),
    current_log_likelihood(0),
    seed(0),
    chunkStepNumber(0),
    nStep(100),
    verbose(false)
{

    if( mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
        int i = 0;
        while( mkdir( (filename+generalUtils::toString(i)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0 ){
            i++;
        }
        filename = filename+generalUtils::toString(i);
    }

    for(int i = 0;i<nVar;i++) trace[i] = std::vector<float>(chunkSize);

    // construct a trivial random generator engine from a time-based seed:
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
}

int main(int argc, char** argv){
    int c;
    int nStep = 0;
    std::string name = "test";
    std::string maskFile = ""; 
    bool verbose = false;
    float alphaRegularization = 0;
    while((c =  getopt(argc, argv, "n:f:v:a:m:")) != EOF)
        {
            switch (c)
                {
                case 'n':
                    nStep = generalUtils::stringTo<int>(optarg);
                    break;
                case 'f':
                    name = optarg;
                    break;
                case 'v':
                    verbose = true;
                    break;
                case 'a':
                    alphaRegularization = generalUtils::stringTo<float>(optarg);
                    break;
                case 'm':
                    maskFile = optarg;
                    break;
                default:
                    break;
                }
        }

    std::clock_t start = std::clock();

    //    MCMC<RealisticToyModel, ProposalFunction > a(name);

    
    RealDataModel model;
    if(maskFile != "") model.SetMask(maskFile);
    model.setRegularizationFactor(alphaRegularization);

    MCMC<RealDataModel, ProposalFunction > a(name,model);
    a.setVerbose(verbose);
    
    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    std::cout << "Processing of : " << name << " done !" << std::endl;
    return 0;
}

