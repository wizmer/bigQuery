#ifndef _MODELBASE__HPP_
#define _MODELBASE__HPP_

template<typename SearchSpaceType>
class ModelBase
{
protected:
    std::vector<float*> spectatorVariables;
    SearchSpaceType realValues;
    SearchSpaceType initialConditions;
public:
    ModelBase() : spectatorVariables(),realValues(), initialConditions(){}

    // ModelBase(const ModelBase & original):
    //     realValues(original.realValues),
    //     initialConditions(original.initialConditions){}

    virtual ~ModelBase(){}

    virtual float GetLogLikelihood(const SearchSpaceType & point) = 0;
    virtual void saveMetaData(const std::ofstream & ){}

    const SearchSpaceType GetRealValues(){ return realValues; }
    const SearchSpaceType GetInitialConditions(){ return initialConditions; }
};

#endif

