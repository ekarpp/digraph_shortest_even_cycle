#ifndef FMATRIX_TEST_H
#define FMATRIX_TEST_H

#include "../src/gf.hh"
#include "../src/fmatrix.hh"
#include "test.hh"

class Matrix_test : Test
{
private:
    int dim;

    void test_addition();
    void test_multiplication();

    FMatrix random_matrix();

    void run()
    {
        test_addition();
    }

public:
    Matrix_test(int dim);
};

#endif
