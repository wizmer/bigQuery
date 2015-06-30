#include <vector>
#include <functional>
#include <iostream>
#include <fstream>

class Matrix 
{
    int nRows,nColumns;
    std::vector<float> data;
public:
    Matrix(int iN, int iM): 
        nRows(iN), nColumns(iM),data(nRows*nColumns,0){}

    inline int getNrows()   const { return nRows;    }
    inline int getNcolums() const { return nColumns; }

    inline float &  at(int n, int m) { return data[n + nRows*m]; }
    inline float   get(int n, int m) const { return data[n + nRows*m]; }
    inline void    set(int n, int m, float val) { data[n + nRows*m] = val; }

    template<typename T>
    void Fill(T matrix)
    {
        for(int n = 0; n < nRows; n++) 
            for(int m = 0; m < nColumns; m++) 
                at(n,m) = matrix[n][m];
    }

    Matrix Dot(const Matrix & rhs) const
    {
        Matrix ret(nRows, rhs.nColumns);
        for(int m = 0; m < rhs.nColumns; m++){
            for(int k = 0; k < nColumns; k++){
                float val = rhs.get(k,m);
                for(int n = 0; n < nRows; n++){
                    ret.at(n,m) += get(n,k)* val;
                }
            }
        }
        return ret;
    }

    Matrix Transpose() const
    {
        Matrix ret(nColumns, nRows);
        for(int m = 0; m < nColumns; m++)
            for(int n = 0; n < nRows; n++) 
                ret.at(m,n) = get(n,m);
        return ret;
    }

    void map(std::function<float(float,int,int)> func)
    {
        for(int m = 0; m < nColumns; m++)
            for(int n = 0; n < nRows; n++) 
                at(n,m) = func(get(n,m), n, m);
    }

    float applyAndSum(std::function<float(float,int,int)> func)
    {
        float ret = 0;
        for(int n = 0; n < nRows; n++) {
            for(int m = 0; m < nColumns; m++) {
                ret += func(get(n,m), n, m);
            }
        }
        return ret;
    }

    void dump()
    {
        for(int n = 0; n < nRows; n++){
            for(int m = 0; m < nColumns; m++){
                std::cout << get(n,m) << "\t";
            }
            std::cout << std::endl;
        }
    }

  void save( std::string filename ){
      std::cout << "filename : " << filename << std::endl;
      std::ofstream f( filename );
      for(int n = 0; n < nRows; n++){
          for(int m = 0; m < nColumns; m++){
              f << get(n,m) << ",";
          }
          f << std::endl;
      }
      f.close();
  }
};

