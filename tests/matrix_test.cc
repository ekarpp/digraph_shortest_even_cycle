#include <iostream>
#include <vector>

#include "matrix_test.hh"
#include "../src/fmatrix.hh"
#include "../src/global.hh"
#include "../src/gf.hh"

using namespace std;

Matrix_test::Matrix_test(int dim)
{
    this->dim = dim;

    cout << "--------------" << endl;
    cout << "TESTING MATRIX" << endl;
    cout << "--------------" << endl << endl;

    this->run();
}

FMatrix Matrix_test::random_matrix()
{
    vector<vector<GF_element>> m(this->dim, vector<GF_element>(this->dim));
    for (int i = 0; i < this->dim; i++)
        for (int j = 0; j < this->dim; j++)
            m[i][j] = global::F.random();
    return FMatrix(this->dim, m);
}

void Matrix_test::test_addition()
{
    cout << "addition" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        FMatrix sum = this->random_matrix();
        for (int i = 0; i < this->dim; i++)
            for (int j = 0; j < this->dim; j++)
                sum.set(i, j, a[i][j] + b[i][j]);

        if (sum != a + b)
            err++;
    }
    end_test(err);
}

void Matrix_test::test_multiplication()
{

}
