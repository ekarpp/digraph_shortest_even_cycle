/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef FMATRIX_TEST_H
#define FMATRIX_TEST_H

#include "test.hh"
#include "../../src/fmatrix.hh"

class FMatrix_test : Test
{
private:
    int dim;

    void test_determinant();

    void test_det_singular();

    void run()
    {
        test_determinant();
        test_det_singular();
    }

    FMatrix vandermonde();
    FMatrix random();

public:
    FMatrix_test(int dim, int tests);
};

#endif
