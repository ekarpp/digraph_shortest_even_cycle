#ifndef GF_TEST_H
#define GF_TEST_H

#include "../src/gf.hh"

class GF_test
{
private:
    GF2n field;
    int n;
    void end_test(int err);
public:
    GF_test(int deg);
    void test_add_inverse();
    void test_mul_commutative();
    void test_mul_id();
    void test_mul_inverse();
};

#endif
