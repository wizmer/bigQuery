
template<typename SearchSpace>
class ModelBase
{
protected:
    std::vector<float*> spectatorVariables;
    SearchSpace realValues;
    SearchSpace initialConditions;
public:
    ModelBase() : spectatorVariables(),realValues(), initialConditions(){}

    // ModelBase(const ModelBase & original):
    //     realValues(original.realValues),
    //     initialConditions(original.initialConditions){}

    virtual ~ModelBase(){}

    virtual float GetLogLikelihood(const SearchSpace & point) = 0;
    virtual void saveMetaData(const std::ofstream & ){}

    const SearchSpace GetRealValues(){ return realValues; }
    const SearchSpace GetInitialConditions(){ return initialConditions; }
};

