#include <vector>
#include <functional>
#include <iostream>
#include <fstream>

class Matrix 
{
    int nRows,nColumns;
    std::vector<float> data;
public:
    Matrix():nRows(0), nColumns(0) {}
    Matrix(int iN, int iM): 
        nRows(iN), nColumns(iM),data(nRows*nColumns,0){}

    inline int getNrows()   const { return nRows;    }
    inline int getNcolums() const { return nColumns; }

    inline float &  at(int row, int column) { return data[row + nRows*column]; }
    inline float   get(int row, int column) const { return data[row + nRows*column]; }
    inline void    set(int row, int column, float val) { data[row + nRows*column] = val; }

    template<typename T>
    void Fill(T matrix)
    {
        for(int iRow = 0; iRow < nRows; iRow++) 
            for(int iColumn = 0; iColumn < nColumns; iColumn++) 
                at(iRow,iColumn) = matrix[iRow][iColumn];
    }

    Matrix Dot(const Matrix & rhs) const
    {
        //if(nColumns != rhs.nRows)
        //{
        //    std::cout << "LHS (" <<     nRows << "," <<     nColumns << ")\n";
        //    std::cout << "RHS (" << rhs.nRows << "," << rhs.nColumns << ")\n";
        //    std::cout << "Dimensions mismatch!\n";
        //    exit(1);
        //}

        Matrix ret(nRows, rhs.nColumns);
        for(int iRhsColumn = 0; iRhsColumn < rhs.nColumns; iRhsColumn++){
            for(int iLhsColumn = 0; iLhsColumn < nColumns; iLhsColumn++){
                float val = rhs.get(iLhsColumn,iRhsColumn);
                for(int n = 0; n < nRows; n++){
                    ret.at(n,iRhsColumn) += get(n,iLhsColumn)* val;
                }
            }
        }
        return ret;
    }

    Matrix Transpose() const
    {
        Matrix ret(nColumns, nRows);
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                ret.at(iColumn,iRow) = get(iRow,iColumn);
        return ret;
    }

    void map(std::function<float(float,int,int)> func)
    {
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                at(iRow,iColumn) = func(get(iRow,iColumn), iRow, iColumn);
    }

    float applyAndSum(std::function<float(float,int,int)> func)
    {
        float ret = 0;
        for(int iRow = 0; iRow < nRows; iRow++) {
            for(int iColumn = 0; iColumn < nColumns; iColumn++) {
                ret += func(get(iRow,iColumn), iRow, iColumn);
            }
        }
        return ret;
    }

    void dump()
    {
        for(int iRow = 0; iRow < nRows; iRow++){
            for(int iColumn = 0; iColumn < nColumns; iColumn++){
                std::cout << get(iRow,iColumn) << "\t";
            }
            std::cout << std::endl;
        }
    }

  void save( std::string filename ){
      std::cout << "filename : " << filename << std::endl;
      std::ofstream f( filename );
      for(int iRow = 0; iRow < nRows; iRow++){
          for(int iColumn = 0; iColumn < nColumns; iColumn++){
              f << get(iRow,iColumn) << ",";
          }
          f << std::endl;
      }
      f.close();
  }
};

