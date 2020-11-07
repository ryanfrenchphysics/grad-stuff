#include <iostream>
#include <functional>
// Use eigen for lin alg
#include <./eigen/Eigen/Dense>
#include <./eigen/Eigen/Eigenvalues>
#include <./eigen/Eigen/Core>
#include <complex>
#include <vector>
#include <iomanip> // Set precisions

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

// Need to use for 'd' and 'id' literals
// (double, imaginary double, respectively)
using namespace std::literals;

void spacer()
{
        std::cout << "-------------------------\n\n";
}

// Test function to map to each element of array
std::complex<double> func(std::complex<double> x)
{
        return conj(x) * x;
}

/*
 *  If we are using std::literals, there's no real need to
 *  use complex, as we see in constructing the matrices
 */
int main()
{
        std::cout << std::setprecision(2);
        std::cout << std::fixed;
        std::cout << std::showpos;
        std::cout <<
        "EIGEN EXAMPLES:\n" <<
        "---------------\n\n";

        // Create 2x2 double matrix
        MatrixXd m1(2, 2);
        // Add values element-wise
        m1(0,0) = 3;
        m1(1,0) = 2.5;
        m1(0,1) = -1;
        m1(1,1) = m1(1,0) + m1(0,1);
        std::cout << "m1:\n" << m1 << "\n\n";
        spacer();

        // Create 3x3 double matrix
        MatrixXd m2(3, 3);
        // Add values with stream operator
        m2 <<   1.0, 2.0, 4,
                2.2, 1.8, 7,
                8.8, 2, 4.3;
        std::cout << "m2:\n" << m2 << "\n\n";
        spacer();


        // Create 2x2 complex double matrix
        MatrixXcd m3(2, 2);
        m3 <<   1.0d + 3.2id, 4.0d + 2.0id,
                3.1d - 3.0id, 2.6d - 2.9id;
        std::cout << "m3:\n" << m3 << "\n\n";
        spacer();

        // apply func to each element in m3
        // Use unaryExpr to map func to each element
        std::cout << "Mod squared each element of m3:\n" << m3.unaryExpr(&func) << "\n\n";
        spacer();

        // Testing complex * conjugate
        std::complex<double> a = (1.0d + 2.0id) * (1.0d - 2.0id);
        std::cout << "(1+2i) * (1-2i): " << a << "\n\n";
        spacer();


        // Print real, imaginary, and complex conjugate values
        std::cout << "Real m3:\n" << m3.real() << "\n" << std::endl;
        std::cout << "Imaginary m3:\n" << m3.imag() << "\n" << std::endl;
        std::cout << "Conj m3:\n" << m3.conjugate() << "\n" << std::endl;
        std::cout << "Adjoint m3:\n" << m3.adjoint() << "\n" << std::endl;
        spacer();

        // Create some complex vectors, size 3
        VectorXcd v1(3);
        v1 <<   1.2d - 3.3id,
                3.2d + 2.0id,
                1.0d - 3.2id;
        std::cout << "v1:\n" << v1 << "\n\n";
        spacer();

        VectorXcd v2(3);
        v2 <<   2.3d + 8.4id,
                1.0d + 2.9id,
                1.1d - 6.7id;

        std::cout << "v1 dot v2:\n" << v1.dot(v2) << "\n\n";
        spacer();

        std::cout << "v1 * m2:\n" << v1.transpose() * m2 << "\n\n";
        spacer();


        // Create object for eigenvals/eigenvecs
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> cev;
        // Get eigenvals, eigenvecs of m2
        cev.compute(m2);
        std::cout << "Eigenvals of m2:\n" << cev.eigenvalues() << "\n\n";
        std::cout << "Eigenvecs of m2:\n" << cev.eigenvectors() << "\n\n";
        spacer();


        // Create the pauli matrices
        MatrixXcd pauli1(2,2);
        MatrixXcd pauli2(2,2);
        MatrixXcd pauli3(2,2);
        MatrixXcd pauli0(2,2);
        pauli0 <<       1.0d + 0.0id, 0.0d + 0.0id,
                        0.0d + 0.0id, 1.0d + 0.0id;

        pauli1 <<       0.0d + 0.0id, 1.0d + 0.0id,
                        1.0d + 0.0id, 0.0d + 0.0id;

        pauli2 <<       0.0d + 0.0id, 0.0d - 1.0id,
                        0.0d + 1.0id, 0.0d + 0.0id;

        pauli3 <<       1.0d + 0.0id, 0.0d + 0.0id,
                        0.0d + 0.0id, -1.0d + 0.0id;

        // Create vector of these matrices
        std::vector<MatrixXcd> pauli_vec(4);
        pauli_vec[0] = pauli0;
        pauli_vec[1] = pauli1;
        pauli_vec[2] = pauli2;
        pauli_vec[3] = pauli3;

        for (int i=0; i<3; i++){
                std::cout << "Pauli matrix " << i+1 << ":\n" << pauli_vec[i] << "\n\n";
        }
        spacer();

        std::cout << "Pauli 1 * Pauli 2:\n" << pauli_vec[0]*pauli_vec[1] << "\n\n";
        spacer();

        // Do eigen operations work with vector subscripting? (YES!)
        std::cout << "Pauli 2 conjugate:\n" << pauli_vec[1].conjugate() << "\n\n";
        spacer();

        return 0;
}

// compile with:
// gcc -lstdc++ -lm -std=c++14 -Wall -fconcepts -fext-numeric-literals -g test.cpp -o test -O[1, 2, 3, s]

// If you have g++, replace "gcc -lstdc++" with g++

/*
 *  timing for optimizations:
 *  None:       18.9s
 *  O1:         24.8s
 *  O2:         29.5s
 *  O3:         31.1s
 *  Os:         23.7s
 */
