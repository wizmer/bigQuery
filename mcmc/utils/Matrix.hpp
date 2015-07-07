#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <functional>
#include <iostream>
#include <fstream>
#include <cmath>

template <typename T> class Matrix 
{
    long unsigned int nRows,nColumns;
    std::vector<T> data;
public:

    Matrix():nRows(0), nColumns(0), data(0) {}
    Matrix(long unsigned int iN, long unsigned int iM): 
        nRows(iN), nColumns(iM),data(nRows*nColumns,0){}

    inline long unsigned int getNrows()   const { return nRows;    }
    inline long unsigned int getNcolums() const { return nColumns; }

    inline T &  at(long unsigned int row, long unsigned int column) { return data[row + nRows*column]; }
    inline T   get(long unsigned int row, long unsigned int column) const { return data[row + nRows*column]; }
    inline void    set(long unsigned int row, long unsigned int column, T val) { data[row + nRows*column] = val; }

    template<typename M>
    void Fill(M matrix)
    {
        for(int iRow = 0; iRow < nRows; iRow++) 
            for(int iColumn = 0; iColumn < nColumns; iColumn++) 
                set(iRow,iColumn, matrix[iRow][iColumn]);
    }

    Matrix<T> Dot(const Matrix & rhs) const
    {
        //if(nColumns != rhs.nRows)
        //{
        //    std::cout << "LHS (" <<     nRows << "," <<     nColumns << ")\n";
        //    std::cout << "RHS (" << rhs.nRows << "," << rhs.nColumns << ")\n";
        //    std::cout << "Dimensions mismatch!\n";
        //    exit(1);
        //}

        Matrix<T> ret(nRows, rhs.nColumns);
        for(int iRhsColumn = 0; iRhsColumn < rhs.nColumns; iRhsColumn++){
            for(int iLhsColumn = 0; iLhsColumn < nColumns; iLhsColumn++){
                T val = rhs.get(iLhsColumn,iRhsColumn);
                for(int n = 0; n < nRows; n++){
                    ret.at(n,iRhsColumn) += get(n,iLhsColumn)* val;
                }
            }
        }
        return ret;
    }

    Matrix<T> Transpose() const
    {
        Matrix<T> ret(nColumns, nRows);
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                ret.at(iColumn,iRow) = get(iRow,iColumn);
        return ret;
    }

    void map(std::function<T(T,int,int)> func)
    {
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                at(iRow,iColumn) = func(get(iRow,iColumn), iRow, iColumn);
    }

    void map(std::function<T(T,int)> func)
    {
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                at(iRow,iColumn) = func(get(iRow,iColumn), iRow);
    }

    void map(std::function<float(float)> func)
    {
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++) 
                at(iRow,iColumn) = func(get(iRow,iColumn));
    }

    void zeroes(){
        std::fill(data.begin(), data.end(), 0);
    }

    void linearCombination(std::vector<Matrix<T>> &matrixBase, std::vector<float> &lambda)
    {
        long unsigned int nMatrices = matrixBase.size();
        for(int iColumn = 0; iColumn < nColumns; iColumn++)
            for(int iRow = 0; iRow < nRows; iRow++)
                for(int iMatrix = 0;iMatrix<nMatrices;iMatrix++)
                    at(iRow,iColumn) += matrixBase[iMatrix].get(iRow,iColumn) * lambda[iMatrix];
    }

    T applyAndSum(std::function<T(T,int,int)> func, const Matrix<bool> & masked )
    {
        T ret = 0;
        for(int iRow = 0; iRow < nRows; iRow++) {
            for(int iColumn = 0; iColumn < nColumns; iColumn++) {
                if( masked.get(iRow, iColumn) == true ) continue;
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

    Matrix<T> subMatrix(long unsigned int _nRows, long unsigned int _nColumns = 0,
                        long unsigned int _firstRow = 0, long unsigned int _firstColumn = 0){
        if( _nRows <= 0 ) _nRows = this -> nRows;
        if( _nColumns <= 0 ) _nColumns = this -> nColumns;

        if( _nRows + _firstRow > nRows || _nColumns + _firstColumn > nColumns ){
            std::cout << "Cannot make a submatrix with such dimensions" << std::endl;
            exit(-1);
        }
        Matrix<T> output(_nRows, _nColumns);
        for(int iRow = 0;iRow<_nRows;iRow++){
            for(int iColumn = 0;iColumn<_nColumns;iColumn++){
                output.set(iRow,iColumn, this -> get(_firstRow + iRow, _firstColumn + iColumn));
            }
	}
        return output;
    }

    inline Matrix<T> operator*(float scale){
        Matrix<T> output(nRows,nColumns);
        for(int iColumn = 0; iColumn < nColumns; iColumn++){
            for(int iRow = 0; iRow < nRows; iRow++){
                output.at(iRow,iColumn) = scale * get(iRow,iColumn);
            }
        }
        return output;
    }


    inline Matrix<T> operator+=(const Matrix &rhs){
        for(int iColumn = 0; iColumn < nColumns; iColumn++){
            for(int iRow = 0; iRow < nRows; iRow++){
                at(iRow,iColumn) += rhs.get(iRow,iColumn);
            }
        }
        return (*this);
    }

    inline bool operator==(const Matrix &rhs){
        if( nRows != rhs.nRows ) return false;
        if( nColumns != rhs.nColumns ) return false;

        for(int iColumn = 0; iColumn < nColumns; iColumn++){
            for(int iRow = 0; iRow < nRows; iRow++){
                if( std::abs(rhs.get(iRow,iColumn) - get(iRow,iColumn)) > 1e-09 ) return false;
            }
        }
        return true;
    }

    inline bool operator!=(const Matrix<T> &rhs){
        return !operator==(rhs);
    }
};

typedef Matrix<float> MatrixF;
typedef Matrix<bool> MatrixB;
#endif //MATRIX_H
