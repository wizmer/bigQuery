#include <iostream>
#include "Matrix.hpp"

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

bool test_copy()
{
    Matrix A(4,3), R(4,3);
    V2D data = { {1,2,3},
                 {4,5,6},
                 {7,8,9},
                 {1,2,3} };
    A.Fill(data);
    R.Fill(data);
 
    Matrix B(4,3);
    B = A;
    A.map([](float v, int n,int m){return 0;});

    bool pass = true;
    B.map([&pass, &R](float v, int n,int m){pass = pass && (R.get(n,m) == v); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " failed\n";
        B.map([](float v, int n,int m){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;

}

bool test_subMatrix(){
    bool pass = true;
    Matrix A(4,3);
    V2D data = { {1,2,3},
                 {4,5,6},
                 {7,8,9},
                 {1,2,3} };
    A.Fill(data);
    Matrix B = A.subMatrix(3,2,1,1);
    if( B.getNrows() != 3) pass = false;
    if( B.getNcolums() != 2) pass = false;
    if( B.get(0,0) != 5 ) pass = false;
    if( B.get(0,1) != 6 ) pass = false;
    if( B.get(1,0) != 8 ) pass = false;
    if( B.get(1,1) != 9 ) pass = false;
    if( B.get(2,0) != 2 ) pass = false;
    if( B.get(2,1) != 3 ) pass = false;

    if(!pass)
        {
            std::cout << __FUNCTION__ << " failed\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " passed.\n";

    return pass;

}
int main(void)
{
    test_dot_1();
    test_dot_2();
    test_copy();
    test_subMatrix();
    return 0;
}
