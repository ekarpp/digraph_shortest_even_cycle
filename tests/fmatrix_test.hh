#ifndef FMATRIX_TEST_H
#define FMATRIX_TEST_H

#include "test.hh"
#include "../src/fmatrix.hh"

class FMatrix_test : Test
{
private:
    int dim;
    int tests = 1;

    void test_determinant();

    void run()
    {
        test_determinant();
    }

    FMatrix vandermonde();

public:
    FMatrix_test(int dim);
};

#endif
