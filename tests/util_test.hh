/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_TEST_H
#define UTIL_TEST_H

#include "test.hh"

class Util_test : Test
{
private:
    int n;

    void test_interpolation();

    void run()
    {
        test_interpolation();
    }

public:
    Util_test(int n, int tests);
};

#endif
