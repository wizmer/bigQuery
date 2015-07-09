#ifndef _SEARCHSPACE__HPP_
#define _SEARCHSPACE__HPP_

struct SearchSpaceBase
{
    SearchSpaceBase(){}
    virtual ~SearchSpaceBase(){}

    virtual long unsigned int size() const = 0;
    virtual float getRaw(int i) const = 0;
};

// std::ostream& operator<<(std::ostream& os, const SearchSpaceBase& point)
// {
//     for(int i = 0;i<point.size();i++){
//         os << point.getRaw(i) << " ";
//     }
//     return os;
// }

struct SearchSpace : public SearchSpaceBase
{
    SearchSpace(): fluxP(), fluxD(){}
    ~SearchSpace(){}

    std::vector<float> fluxP;
    std::vector<float> fluxD;
    inline long unsigned int size() const {return fluxP.size() + fluxD.size();} 
    inline float getRaw(int i) const { return i >= fluxP.size()? fluxD[i-fluxP.size()]:fluxP[i]; } 
};

struct SearchSpace3Param : public SearchSpaceBase
{
    SearchSpace3Param(): fluxP(), fluxD(), fluxK(){}
    ~SearchSpace3Param(){}

    std::vector<float> fluxP;
    std::vector<float> fluxD;
    std::vector<float> fluxK;
    inline long unsigned int size() const {return fluxP.size() + fluxD.size() + fluxK.size();} 
    inline float getRaw(int i) const { 
        if(i < fluxP.size()) return fluxP[i];
        else if ( i < fluxP.size() + fluxD.size() ) return fluxD[i-fluxP.size()];
        else return fluxD[i-fluxP.size()-fluxD.size()];
    } 
};

#endif

