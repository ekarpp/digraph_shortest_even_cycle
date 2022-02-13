#ifndef EXTENSION_TEST_H
#define EXTENSION_TEST_H

#include "../src/extension.hh"

class Extension_test
{
private:
    int n;
    int tests = 10000;
    void end_test(int err);

public:
    Extension_test(int deg);
    void test_add_inverse();
    void test_associativity();
    void test_mul();
};

#endif
