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
        for (int x = 0; x < this->dim; x++)
            for (int y = 0; y < this->dim; y++)
                sum.set(x, y, a[x][y] + b[x][y]);

        if (sum != a + b)
            err++;
    }
    end_test(err);
}

void Matrix_test::test_subtraction()
{
    cout << "subtraction" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        FMatrix diff = this->random_matrix();
        for (int x = 0; x < this->dim; x++)
            for (int y = 0; y < this->dim; y++)
                diff.set(x, y, a[x][y] - b[x][y]);

        if (diff != a - b)
            err++;
    }
    end_test(err);
}

void Matrix_test::test_multiplication()
{
    cout << "multiplication" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        FMatrix prod = this->random_matrix();
        for (int x = 0; x < this->dim; x++)
        {
            for (int y = 0; y < this->dim; y++)
            {
                GF_element sum = global::F.zero();
                for (int j = 0; j < this->dim; j++)
                    sum += a[x][j] * b[j][y];
                prod.set(x, y, sum);
            }
        }

        if (prod != a*b)
            err++;
    }
    end_test(err);
}
