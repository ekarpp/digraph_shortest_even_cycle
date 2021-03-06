/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_TEST_H
#define EXTENSION_TEST_H

#include "test.hh"

class Extension_test : Test
{
private:
    int n;

    void test_add_inverse();
    void test_associativity();
    void test_mul();
    void test_fast_mul();
    void test_intel_rem();
    void test_mont_rem();
    void test_even_tau();
    void test_is_even();
    void test_kronecker_mul();

#if GF2_bits == 16
    void test_packed_intel_rem();
    void test_packed_fast_mul();
#endif

    void run()
    {
        test_add_inverse();
        test_associativity();
        test_mul();
        test_even_tau();
        test_is_even();
        test_fast_mul();
        test_intel_rem();
        test_mont_rem();
        test_kronecker_mul();
#if GF2_bits == 16
        test_packed_intel_rem();
        test_packed_fast_mul();
#endif
    }

public:
    Extension_test(int tests);

};

#endif
