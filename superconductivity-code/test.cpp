// #include "Grid.hpp"
// #include "SCUniversal.hpp"
#include <iostream>
// Use eigen for lin alg
#include <./eigen/Eigen/Dense>
#include <complex>

using Eigen::MatrixXd;
using Eigen::MatrixXcd;

// Need to use for 'd' and 'id' literals
// (double, imaginary double, respectively)
using namespace std::literals;

/*
 *  If we are using std::literals, there's no real need to
 *  use complex, as we see in constructing the matrices
 */
int main()
{
        // Create 2x2 double matrix
        MatrixXd m(2, 2);
        // Add values element-wise
        m(0,0) = 3;
        m(1,0) = 2.5;
        m(0,1) = -1;
        m(1,1) = m(1,0) + m(0,1);
        std::cout << m << "\n" << std::endl;

        // Create 3x3 double matrix
        MatrixXd m2(3, 3);
        // Add values with stream operator
        m2 <<   1.0, 2.0, 4,
                2.2, 1.8, 7,
                8.8, 2, 4.3;
        std::cout << m2 << "\n" << std::endl;

        // Create 2x2 complex double matrix
        MatrixXcd m3(2, 2);
        m3 <<   1.0d + 3.2id, 4.0d + 2.0id,
                3.1d - 3.0id, 2.6d - 2.9id;
        std::cout << m3 << "\n" << std::endl;


        // Print real, imaginary, and complex conjugate values
        std::cout << "Real:\n" << m3.real() << "\n" << std::endl;
        std::cout << "Imaginary:\n" << m3.imag() << "\n" << std::endl;
        std::cout << "Conj:\n" << m3.conjugate() << "\n" << std::endl;


        return 0;
}

// compile with:
// gcc -lstdc++ -lm -std=c++14 -Wall -fconcepts -fext-numeric-literals -g eigen_ex.cpp -o eigen_ex
