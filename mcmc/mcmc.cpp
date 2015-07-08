#include "mcmc.hpp"

template <class ProposalFunction > MCMC<ProposalFunction >::MCMC(std::string name, ModelBase & _model):
    model(_model),
    realValues(model.GetRealValues()),
    current_point(model.GetInitialConditions()),
    initialConditions(model.GetInitialConditions()),
    proposalFunction(model.GetInitialConditions()),
    nVar(initialConditions.size()),
    nThreads(1),
    chunkSize(maxRAM / ( (nVar+1)*sizeof(float)) / nThreads),
    spectators(),
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

template <class ProposalFunction > void MCMC<ProposalFunction >::loop(){
    std::cout << "nStep : " << nStep/1e6 << " millions" << std::endl;
    std::cout << "chunkSize : " << chunkSize/1e6 << " MBytes" << std::endl;
    std::cout << "maxRAM : " << maxRAM/1e6 << " MBytes" << std::endl;
    std::cout << "nVar : " << nVar << std::endl;
    std::cout << "sizeof(float) : " << sizeof(float) << std::endl;

    sleep(1);

    saveMetaData();

    current_log_likelihood = model.GetLogLikelihood(initialConditions);
        
    for( int i = 0; i < nStep; i++)
        {
            if( i%30000 == 0) printf("%i/%i : %i%%\n",i, nStep, int(float(i)/float(nStep)*100));

            auto proposed_point = proposalFunction.proposePoint(current_point);

            float proposed_log_likelihood = model.GetLogLikelihood(proposed_point);

            float the_likelihood_ratio = exp(proposed_log_likelihood-current_log_likelihood);

            if( verbose ){

                std::cout << "fluxP "; for(auto v: proposed_point.fluxP) std::cout << v << " "; std::cout << "\n";
                std::cout << "fluxD "; for(auto v: proposed_point.fluxD) std::cout << v << " "; std::cout << "\n";

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

template <class ProposalFunction > void MCMC<ProposalFunction>::setSpectator(std::function<float(const ModelBase&)> f){
    spectators.push_back(f);
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

    MCMC<ProposalFunction > a(name,model);
    a.setVerbose(verbose);
    //    a.setSpectator( &RealDataModel::regularizationTermAccessor );

    if( nStep > 0 ) a.setSteps(nStep);
    a.run();
    std::cout << "sizeof(a) : " << sizeof(a) << std::endl;
    std::cout << "Time : " << (std::clock() - start) / (float)(CLOCKS_PER_SEC) << " s" << std::endl;
    std::cout << "Processing of : " << name << " done !" << std::endl;
    return 0;
}

