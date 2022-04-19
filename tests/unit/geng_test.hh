/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GENG_TEST_H
#define GENG_TEST_H

#include "test.hh"

class Geng_test : Test
{
private:
    int n;

    void test_geng();

    void run()
    {
        test_geng();
    }

public:
    Geng_test();
};

#endif
