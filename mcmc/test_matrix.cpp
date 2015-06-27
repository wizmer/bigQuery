#include <iostream>
#include "pd_model.hpp"

typedef std::vector<std::vector<float>> V2D;

bool test_dot_1()
{
    Matrix A(4,3);
    V2D dataA = { {1,2,3},
                  {4,5,6},
                  {7,8,9},
                  {1,2,3} };
    A.Fill(dataA);

    Matrix B(3,3);
    V2D dataB = { {1,0,0},
                  {0,1,0},
                  {0,0,1} };
    B.Fill(dataB);

    Matrix C = A.Dot(B);

    bool pass = true;
    C.map([&pass, &A](float v, int n,int m){pass = pass && (A.get(n,m) == v); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " failed:\n";
        C.map([](float v, int n,int m){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;
}

bool test_dot_2()
{
    Matrix A(2,2);
    V2D dataA = {{2,4}, {3,1}};
    A.Fill(dataA);

    Matrix B(2,2);
    V2D dataB = {{3,1}, {2,-1}};
    B.Fill(dataB);

    Matrix R(2,2);
    V2D dataR = { { 2*3 + 4*2, 2*1 + 4*(-1)}, 
                  { 3*3 + 1*2, 3*1 + 1*(-1)}};
    R.Fill(dataR);

    Matrix C = A.Dot(B);

    bool pass = true;
    C.map([&pass, &R](float v, int n,int m){pass = pass && (R.get(n,m) == v); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " failed\n";
        C.map([](float v, int n,int m){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;
}
int main(void)
{
    test_dot_1();
    test_dot_2();
    return 0;
}
