/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <valarray>

#include "ematrix_test.hh"
#include "../../src/global.hh"
#include "../../src/ematrix.hh"
#include "../../src/extension.hh"

using namespace std;

EMatrix_test::EMatrix_test(int tests = 0)
{
    if (tests)
        this->tests = tests;

    cout << "---------------" << endl;
    cout << "TESTING EMATRIX" << endl;
    cout << "---------------" << endl;

    this->run();
}

EMatrix EMatrix_test::random()
{
    valarray<Extension_element> m(this->dim * this->dim);

    for (int row = 0; row < this->dim; row++)
        for (int col = 0; col < this->dim; col++)
            m[row*this->dim + col] = global::E.random();

    return EMatrix(this->dim, m);
}

Extension_element EMatrix_test::term(valarray<int> &perm, EMatrix &m)
{
    Extension_element ret = global::E.one();
    for (int col = 0; col < this->dim; col++)
        ret *= m(perm[col], col);
    return ret;
}

void EMatrix_test::swap(int i1, int i2, valarray<int> &perm)
{
    int tmp = perm[i1];
    perm[i1] = perm[i2];
    perm[i2] = tmp;
}

void EMatrix_test::test_per_det()
{
    cout << "per minus det: ";
    int err = 0;
    for (int t = 0; t < this->tests; t++)
    {
        EMatrix m = this->random();
        Extension_element per = global::E.zero();
        Extension_element det = global::E.zero();

        /* iterative heaps algo for permutations. compute permanent
         * and determinant with the Leibniz formula */
        valarray<int> c(0, this->dim);
        valarray<int> perm(0, this->dim);
        for (int i = 0; i < this->dim; i++)
            perm[i] = i;
        Extension_element tt = this->term(perm, m);
        per += tt;
        det += tt;

        bool neg = true;
        int i = 0;
        while (i < this->dim)
        {
            if (c[i] < i)
            {
                if (i%2 == 0)
                    this->swap(0, i, perm);
                else
                    this->swap(c[i], i, perm);
                tt = this->term(perm, m);
                per += tt;
                if (neg)
                    det -= tt;
                else
                    det += tt;
                neg = !neg;
                c[i]++;
                i = 0;
            }
            else
            {
                c[i] = 0;
                i++;
            }
        }
        /* heaps algo end */

        if (per - det != m.per_m_det())
            err++;
    }
    end_test(err);
}
