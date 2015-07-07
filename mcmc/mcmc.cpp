#include "mcmc.hpp"


template <class Model, class ProposalFunction > MCMC<Model, ProposalFunction >::MCMC(std::string name):
    model(),
    realValues(model.realValues),
    current_point(model.initialConditions),
    initialConditions(model.initialConditions),
    proposalFunction(model.initialConditions),
    nVar(initialConditions.size()),
    log_likelihood(chunkSize),
    trace(nVar),
    filename(name),
    current_log_likelihood(0),
    regularizationFactor(0),
    seed(0),
    chunkStepNumber(0),
    nThreads(1), nStep(100),
    chunkSize(maxRAM / ( (nVar+1)*sizeof(float)) / nThreads),
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


template <class Model, class ProposalFunction > float MCMC<Model, ProposalFunction >::GetLogLikelihood(const SearchSpace & point){
    float proposed_log_likelihood = model.GetLogLikelihood(point);
    float antismoothness = model.GetRegularizationTerm(point);
    // std::cout << "regularizationFactor : " << regularizationFactor << std::endl;
    // std::cout << proposed_log_likelihood << "\t" << antismoothness*regularizationFactor << std::endl;
    proposed_log_likelihood -= regularizationFactor * antismoothness;
    return proposed_log_likelihood;
}

int main(int argc, char** argv){
    int c;
    int nStep = 0;
    std::string name = "test";
    bool verbose = false;
    float alphaRegularization = 0;

    while((c =  getopt(argc, argv, "n:f:v:a:")) != EOF)
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
                default:
                    break;
                }
        }

    std::clock_t start = std::clock();

    //    MCMC<RealisticToyModel, ProposalFunction > a(name);
    MCMC<RealDataModel, ProposalFunction > a(name);
    a.setRegularizationFactor(alphaRegularization);
    a.setVerbose(verbose);
    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    std::cout << "Processing of : " << name << " done !" << std::endl;
    return 0;
}

