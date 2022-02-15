/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <valarray>

#include "fmatrix_test.hh"
#include "../src/global.hh"
#include "../src/fmatrix.hh"
#include "../src/gf.hh"

using namespace std;

FMatrix_test::FMatrix_test(int dim, int tests = 0)
{
    if (tests)
        this->tests = tests;
    this->dim = dim;

    cout << "---------------" << endl;
    cout << "TESTING FMATRIX" << endl;
    cout << "---------------" << endl;

    this->run();
}
#include <bitset>
FMatrix FMatrix_test::vandermonde()
{
    valarray<GF_element> m(this->dim * this->dim);
    uint64_t v = global::F.rem(global::randgen());

    for (int row = 0; row < this->dim; row++)
    {
        const GF_element e = GF_element(v);
        v = global::F.rem(v + 2);
        GF_element prod = global::F.one();

        for (int col = 0; col < this->dim; col++)
        {
            m[row*this->dim + col] = prod;
            prod *= e;
        }
    }

    return FMatrix(this->dim, m);
}

void FMatrix_test::test_determinant()
{
    cout << "determinant" << endl;
    int err = 0;
    for (int t = 0; t < this->tests; t++)
    {
        FMatrix vander = this->vandermonde();
        GF_element det = global::F.one();

        for (int i = 0; i < this->dim; i++)
            for (int j = i+1; j < this->dim; j++)
                det *= vander(j, 1) - vander(i, 1);

        if (det != vander.det())
            err++;
    }
    end_test(err);
}
