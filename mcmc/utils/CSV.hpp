#include "generalUtils.hpp"

std::vector<float> getRow(float & first, std::string str)
{
    std::stringstream ss(str);
    std::vector<float> row;
    std::string entry;
    float value;

    getline(ss, entry, ',');

    {
        std::stringstream se(entry);
        if(!(se >> first)) first = 0;
    }

    while (getline(ss, entry, ','))
    {
        std::stringstream se(entry);
        if(!(se >> value)) value = 0;
        row.push_back(value);
    }
    return row;
}

MatrixF getMatrixAndBins( std::fstream & fs, 
                          std::vector<float>& binsT, std::vector<float>& binsM )
{
    float first;
    std::string line;
    std::vector< std::vector<float> > data;

    getline(fs,line);
    binsM = getRow(first, line);

    while (getline(fs,line))
        {
            data.push_back(getRow(first, line));
            binsT.push_back(first);
        }

    MatrixF M(binsT.size()-1, binsM.size()-1);
    M.Fill(data);

    std::vector<float> sums;
    for(auto row : data)
        {
            float sum = 0;
            for(auto v : row) sum += v;
            sums.push_back(sum);
        }

    M.map([&sums](float v, int t ){return v/(sums[t]>0?sums[t]:1);});


    return M;
}

MatrixB getMask(const std::string & fname){
    std::vector<std::string> rows = generalUtils::splitIntoLines( generalUtils::exec("cat " + fname) );
    long unsigned int nRows = rows.size();
    if(nRows == 0){
        std::cerr << "WARNING in CSV.hpp:getMask" << std::endl;
        std::cerr << "Unable to read file: " << fname << std::endl;
        std::cerr << "Exit !" << std::endl;
        exit(-1);
    }

    int nColumns = 0;
    std::vector< std::vector<bool> > data;
    for(int i = 0;i<rows.size();i++){
        std::vector<bool> vec = generalUtils::stringTo<bool>(generalUtils::split(rows[i],","));
        if( i == 0) nColumns = vec.size();
        else if( vec.size() != nColumns ){
            std::cerr << "WARNING in CSV.hpp:getMask" << std::endl;
            std::cerr << "Column mismatch in file: " << fname << std::endl;
            std::cerr << "Exit !" << std::endl;
            exit(-1);
        }
        data.push_back(vec);
    }
    
    MatrixB output(nRows, nColumns);
    output.Fill<std::vector< std::vector<bool> > >(data);
    return output;
}

