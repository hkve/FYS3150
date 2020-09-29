#include <iostream>
#include <cmath>
#include <iomanip>
#include "JacobiEigSolver.hpp"

using namespace std;

void set3by3TestMatrix(double** A, double* eigvals, double** eigvecs);
double dot(double* v1, double* v2, int N);
double offdiag(int N, double** A);
bool TestOrthonorm(int N, double** testMatrix);
bool TestConvergence(int N, double** testMatrix);

int main(int argc, char** argv)
{
    int no_PASSED = 0; // passed tests
    int N; // test matrix dimensions
    double** testMatrix; // the test matrix
    double* ana_eigvals; // analytical eigenvalues of test matrix
    double** ana_eigvecs; // analytical eigenvectors of test matrix

    N = 3;
    ana_eigvals = new double [N];
    testMatrix = new double* [N];
    ana_eigvecs = new double* [N];

    for (int i = 0; i < N; i++) {
        testMatrix[i] = new double [N];
        ana_eigvecs[i] = new double [N];
    }


    // test orthonormality
    set3by3TestMatrix(testMatrix, ana_eigvals, ana_eigvecs);
    bool ORTHONORMAL = TestOrthonorm(N, testMatrix);
    if (ORTHONORMAL) {
        no_PASSED += 1;
    }

    // test convergence
    set3by3TestMatrix(testMatrix, ana_eigvals, ana_eigvecs);
    bool CONVERGE = TestConvergence(N, testMatrix);
    if (CONVERGE) {
        no_PASSED += 1;
    }

    cout << "Tests completed, passed " << no_PASSED << "/2 tests." << endl;
}



void set3by3TestMatrix(double** A, double* eigvals, double** eigvecs) {
    // setting the diagonal elements
    A[0][0] = 3;
    A[1][1] = 6;
    A[2][2] = 3;
    // setting the triangular elements
    A[0][1] = A[1][0] = -2;
    A[0][2] = A[2][0] = 4;
    A[1][2] = A[2][1] = 2;
    
    // setting the analytical solutions
    eigvals[0] = -2;
    eigvals[1] = 7;
    eigvals[2] = 7;

    eigvecs[0][0] = -2.0/3.0; eigvecs[1][0] = -1.0/3.0; eigvecs[2][0] = 2.0/3.0;
    eigvecs[0][1] = -1.0/sqrt(5.0); eigvecs[1][1] = 2.0/sqrt(5.0); eigvecs[2][1] = 0.0;
    eigvecs[0][2] = 1.0/sqrt(2.0); eigvecs[1][2] = 0.0; eigvecs[2][2] = 1.0/sqrt(2.0);
}



bool TestOrthonorm(int N, double** testMatrix) {
    bool PASSED = true;
    double tolerance = 1e-10;
    int k, l;
    double** num_eigvecs;

    JacobiEigSolver* test_solver = new JacobiEigSolver(testMatrix, N, true);
    test_solver -> setTolerance(tolerance);
    test_solver -> Solve();

    num_eigvecs = test_solver -> getEigvecs();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double* eigvec_i = new double [N];
            double* eigvec_j = new double [N];
            for (int k = 0; k < N; k++) {
                eigvec_i[k] = num_eigvecs[k][i];
                eigvec_j[k] = num_eigvecs[k][j];
            }
            double dotprod = dot(eigvec_i, eigvec_j, N);
            if (i == j) {
                if (fabs(dotprod - 1) > tolerance) {
                    cout << "TEST FAILED: The eigenvectors were not normalized!" << endl;
                    cout << "Vector that failed was eigenvector " << i << ":" << endl;
                    cout << "\t";
                    for (int k = 0; k < N; k++) {
                        cout << eigvec_i[k] << ", ";
                    }
                    cout << endl;
                    PASSED = false;
                }
            } else {
                if (fabs(dotprod) > tolerance) {
                    cout << "The eigenvectors were not orthogonal!" << endl;
                    cout << "Vectors that failed were eigenvectors " << i << "and " << j << ":" << endl;
                    cout << "\t";
                    for (int k = 0; k < N; k++) {
                        cout << eigvec_i[k] << ", ";
                    }
                    cout << endl << "\t";
                    for (int k = 0; k < N; k++) {
                        cout << eigvec_i[k] << ", ";
                    }
                    cout << endl;
                    PASSED = false;
                }
            }
        }
    }
    if (PASSED) {
        cout << "TEST PASSED: The eigenvectors were orthonormal." << endl;
    }
    return PASSED;
}



bool TestConvergence(int N, double** testMatrix) {
    bool PASSED = true;
    double max;
    double** resultMatrix;
    double test_offdiag, result_offdiag;

    test_offdiag = offdiag(N, testMatrix);

    resultMatrix = new double* [N];
    for (int i = 0; i < N; i++) {
        resultMatrix[i] = new double [N];
    }
    int k, l;
    JacobiEigSolver* test_solver = new JacobiEigSolver(testMatrix, N, true);

    for (int count = 1; count <= 3; count++) {
        max = 0.0;
        test_solver -> getMax_(&max, &k, &l);
        test_solver -> doJacobiRotation_(k, l);
        
        resultMatrix = test_solver -> getA();
        result_offdiag = offdiag(N, resultMatrix);

        if (result_offdiag > test_offdiag) {
            cout << "TEST FAILED: The algorithms does not converge. The original offdiag was ";
            cout << test_offdiag << " before rotation and " << result_offdiag << " afterwards." << endl;
            PASSED = false;
        }
    }


    if (PASSED) {
        cout << "TEST PASSED: The algorithm converges." << endl;
    }
    return PASSED;
}



double dot(double* v1, double* v2, int N) {
    double product = 0.0;
    for (int i = 0; i < N; i++) {
        product += v1[i] * v2[i];
    }
    return product;
}



double offdiag(int N, double** A) {
    // Compute the square sum of the offdiagonal elements of symmetric NxN matrix A.
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            result += 2 * A[i][j]*A[i][j];
        }
    }
    return sqrt(result);
}