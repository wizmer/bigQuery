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
    C.map([&pass, &A](float v, int n,int m){pass = pass && (std::abs(A.get(n,m) - v) < 1e-99); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m:\n";
        C.map([](float v){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

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
    C.map([&pass, &R](float v, int n,int m){pass = pass && (std::abs(R.get(n,m) - v) < 1e-99); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
        C.map([](float v){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

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
    A.zeroes();

    bool pass = true;
    B.map([&pass, &R](float v, int n,int m){pass = pass && (std::abs(R.get(n,m) - v) < 1e-99); return v;});

    if(!pass)
    {
        std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
        B.map([](float v){std::cout << v << ","; return v;});
        std::cout << "\n";
    }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

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
    if( std::abs(B.get(0,0) - 5) > 1e-99 ) pass = false;
    if( std::abs(B.get(0,1) - 6) > 1e-99 ) pass = false;
    if( std::abs(B.get(1,0) - 8) > 1e-99 ) pass = false;
    if( std::abs(B.get(1,1) - 9) > 1e-99 ) pass = false;
    if( std::abs(B.get(2,0) - 2) > 1e-99 ) pass = false;
    if( std::abs(B.get(2,1) - 3) > 1e-99 ) pass = false;

    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

    return pass;

}

bool test_operatorEqual(){
    V2D data = { {1,2,3},
                 {4,5,6},
                 {7,8,9},
                 {1,2,3} };

    V2D dataWrong = { {1,2,3},
                      {4,7,6},
                      {7,8,9},
                      {1,2,3} };

    Matrix A(4,3),B(4,3),C(4,3);
    A.Fill(data);
    B.Fill(data);
    C.Fill(dataWrong);

    bool pass = true;
    if( (A == B) == false ) pass = false;
    if( A == C ) pass = false;

    if( A != B ) pass = false;
    if( (A != C) == false ) pass = false;

    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

    return pass;
}

bool test_operatorMultiply(){
    V2D data = { {1,2,3},
                 {4,5,6},
                 {7,8,9},
                 {1,2,3} };

    V2D dataTimes2 = { {2,4,6},
                       {8,10,12},
                       {14,16,18},
                       {2,4,6} };

    Matrix A(4,3),B(4,3);
    A.Fill(data);
    B.Fill(dataTimes2);

    bool pass = true;
    Matrix C = A*2;
    if( C != B) pass = false;

    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

    return pass;
}

bool test_operatorAdd(){
    V2D dataA = { {1,2,3},
                 {4,5,6},
                 {7,8,9},
                 {1,2,3} };

    V2D dataB = { {2,4,6},
                  {8,10,12},
                  {14,16,18},
                  {2,4,6} };

    V2D dataC = { {3,6,9},
                  {12,15,18},
                  {21,24,27},
                  {3,6,9} };
    
    Matrix A(4,3),B(4,3), C(4,3);
    A.Fill(dataA);
    B.Fill(dataB);
    C.Fill(dataC);
    A += B;
    
    bool pass = true;
    if( A != C) pass = false;

    if(!pass)
        {
            std::cout << __FUNCTION__ << " \033[1;31mFAILED\033[0m\n";
            std::cout << "\n";
        }
    else std::cout << __FUNCTION__ << " \033[1;32mPASSED\033[0m.\n";

    return pass;
}


int main(void)
{
    test_dot_1();
    test_dot_2();
    test_copy();
    test_subMatrix();
    test_operatorEqual();
    test_operatorMultiply();
    test_operatorAdd();

    return 0;
}
