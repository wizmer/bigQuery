#include <vector>
#include <functional>
#include <iostream>

class Matrix 
{
    int N,M;
    std::vector<float> data;
public:
    Matrix(int iN, int iM): 
        N(iN), M(iM),data(N*M,0){}

    inline float &  at(int n, int m) { return data[n + N*m]; }
    inline float   get(int n, int m) const { return data[n + N*m]; }
    inline void    set(int n, int m, float val) { data[n + N*m] = val; }

    template<typename T>
    void Fill(T matrix)
    {
        for(int n = 0; n < N; n++) 
            for(int m = 0; m < M; m++) 
                at(n,m) = matrix[n][m];
    }

    Matrix Dot(const Matrix & rhs) const
    {
        Matrix ret(N, rhs.M);
        for(int n = 0; n < N; n++){
	  for(int m = 0; m < rhs.M; m++){
	    for(int k = 0; k < M; k++){
	      ret.at(n,m) += get(n,k)*rhs.get(k,m);
	    }
	  }
	}
	return ret;
    }

    void map(std::function<float(float,int,int)> func)
    {
      for(int m = 0; m < M; m++)
        for(int n = 0; n < N; n++) 
	  at(n,m) = func(get(n,m), n, m);
    }

    float applyAndSum(std::function<float(float,int,int)> func)
    {
        float ret = 0;
        for(int n = 0; n < N; n++) {
	  for(int m = 0; m < M; m++) {
	    //	    std::cout << n << "\t" << m << "\t" << get(n,m) << std::endl;
	    ret += func(get(n,m), n, m);
	  }
	}
        return ret;
    }

    void dump()
    {
        for(int n = 0; n < N; n++){
            for(int m = 0; m < M; m++){
                std::cout << get(n,m) << "\t";
            }
            std::cout << std::endl;
        }
    }
};

