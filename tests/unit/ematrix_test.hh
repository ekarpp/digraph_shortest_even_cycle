/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EMATRIX_TEST_H
#define EMATRIX_TEST_H

#include <valarray>

#include "test.hh"
#include "../../src/ematrix.hh"

class EMatrix_test : Test
{
private:
    int dim = 2;

    void test_per_det();

    void run()
    {
        test_per_det();
    }

    EMatrix random();
    Extension_element term(std::valarray<int> &perm, EMatrix &m);
    void swap(int i1, int i2, std::valarray<int> &perm);

public:
    EMatrix_test(int tests);
};

#endif
