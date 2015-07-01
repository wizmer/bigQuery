
std::vector<float> getRow(float & first, std::string str)
{
    std::stringstream ss(str);
    std::vector<float> row;
    std::string entry;
    float value;

    getline(ss, entry, ',');
    std::stringstream se(entry);
    if(!(se >> first)) first = 0;

    while (getline(ss, entry, ','))
    {
        std::stringstream se(entry);
        if(!(se >> value)) value = 0;
        row.push_back(value);
    }
    return row;
}

Matrix getMatrixAndBins( std::fstream & fs, 
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

    Matrix M(binsT.size()-1, binsM.size()-1);
    M.Fill(data);

    std::vector<float> sums;
    for(auto row : data)
    {
        float sum = 0;
        for(auto v : row) sum += v;
        sums.push_back(sum);
    }

    M.map([&sums](float v, int t, int m){return v/(sums[t]>0?sums[t]:1);});

    return M;
}
