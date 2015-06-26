#include <iostream>
#include "pd_model.hpp"


bool test_dot_1()
{

    Matrix A(4,3);
    std::vector<std::vector<float>> dataA = {
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {1,2,3}
    };
    A.Fill(dataA);


    Matrix B(3,3);
    std::vector<std::vector<float>> dataB = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
    B.Fill(dataB);

    Matrix C = A.Dot(B);
    C.map([](int n,int m,float v){std::cout << v << ","; return v; });

    return true;
}


int main(void)
{
    test_dot_1();
    return 0;
}
