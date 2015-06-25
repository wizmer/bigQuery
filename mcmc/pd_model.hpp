#include <vector>
#include <functional>

class Matrix 
{
    int N,M;
    std::vector<float> data;
public:
    Matrix(int iN, int iM): 
        N(iN), M(iM),data(std::vector<float>(N*M,0)){}

    inline float &  at(int n, int m) { return data[n + N*m]; }
    inline float   get(int n, int m) const { return data[n + N*m]; }

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
        for(int n = 0; n < N; n++) 
            for(int m = 0; m < rhs.M; m++) 
                for(int k = 0; k < M; k++) 
                    ret.at(n,m) += get(n,k)*rhs.get(k,m);
    }

    void map(std::function<float(float,int,int)> func)
    {
        for(int n = 0; n < N; n++) 
            for(int m = 0; m < M; m++) 
                at(n,m) = func(get(n,m), n, m);
    }
};


class PDModel
{
    std::vector<float> betaBinsT, betaBinsM;
    std::vector<float> rgdtBinsT, rgdtBinsM;
    Matrix  rgdtF,  betaF;
    Matrix deltaP, deltaD;

public:
    static const float mp;
    static const float md;

    PDModel( const std::vector<float> & bT, const std::vector<float> & bM, 
             const std::vector<float> & rT, const std::vector<float> & rM  );

    Matrix GetPrediction( const std::vector<float> & fluxP,
                          const std::vector<float> & fluxD  );

    template<typename T> 
    void SetRigidityResolution(T matrix){ rgdtF.Fill(matrix); }

    template<typename T>
    void SetBetaResolution    (T matrix){ betaF.Fill(matrix); }
};
