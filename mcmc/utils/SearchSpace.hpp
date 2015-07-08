#ifndef _SEARCHSPACE__HPP_
#define _SEARCHSPACE__HPP_

struct SearchSpace
{
    SearchSpace(): fluxP(), fluxD(){}
    ~SearchSpace(){}

    std::vector<float> fluxP;
    std::vector<float> fluxD;
    inline long unsigned int size(){return fluxP.size() + fluxD.size();}
    inline float getRaw(int i){ return i >= fluxP.size()? fluxD[i-fluxP.size()]:fluxP[i]; }
};

#endif

