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

    void run()
    {
        test_determinant();
    }

    FMatrix vandermonde();

public:
    FMatrix_test(int dim, int tests);
};

#endif
