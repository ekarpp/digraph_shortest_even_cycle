/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <valarray>

#include "matrix_test.hh"
#include "../src/fmatrix.hh"
#include "../src/global.hh"
#include "../src/gf.hh"

using namespace std;

Matrix_test::Matrix_test(int dim, int tests = 0)
{
    this->dim = dim;
    if (tests)
        this->tests = tests;
    cout << "--------------" << endl;
    cout << "TESTING MATRIX" << endl;
    cout << "--------------" << endl;

    this->run();
}

FMatrix Matrix_test::random_matrix()
{
    valarray<GF_element> m(this->dim * this->dim);
    for (int i = 0; i < this->dim; i++)
        for (int j = 0; j < this->dim; j++)
            m[i*this->dim + j] = global::F.random();
    return FMatrix(this->dim, m);
}

void Matrix_test::test_addition()
{
    cout << "addition: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        FMatrix sum = this->random_matrix();
        for (int x = 0; x < this->dim; x++)
            for (int y = 0; y < this->dim; y++)
                sum.set(x, y, a(x,y) + b(x,y));

        if (sum != a + b)
            err++;
    }
    end_test(err);
}

void Matrix_test::test_subtraction()
{
    cout << "subtraction: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        FMatrix diff = this->random_matrix();
        for (int x = 0; x < this->dim; x++)
            for (int y = 0; y < this->dim; y++)
                diff.set(x, y, a(x,y) - b(x,y));

        if (diff != a - b)
            err++;
    }
    end_test(err);
}

void Matrix_test::test_multiplication()
{
    cout << "multiplication: ";
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
                    sum += a(x,j) * b(j,y);
                prod.set(x, y, sum);
            }
        }

        if (prod != a*b)
            err++;
    }
    end_test(err);
}
