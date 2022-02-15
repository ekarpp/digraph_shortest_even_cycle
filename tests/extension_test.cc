/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "extension_test.hh"
#include "../src/extension.hh"
#include "../src/util.hh"
#include "../src/global.hh"

using namespace std;


Extension_test::Extension_test(int deg)
{
    this->n = deg;

    cout << "-----------------" << endl;
    cout << "TESTING EXTENSION" << endl;
    cout << "-----------------" << endl;

    this->run();
}

void Extension_test::test_add_inverse()
{
    cout << "add inverse" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element e = global::E.random();
        if (e - e != global::E.zero())
            err++;
    }
    end_test(err);
}

void Extension_test::test_associativity()
{
    cout << "test associativity" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();
        Extension_element c = global::E.random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    end_test(err);
}

void Extension_test::test_mul()
{
    cout << "test mul" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        if (a*b != b*a || a*global::E.one() != a
            || b*global::E.one() != b
            || a*global::E.zero() != global::E.zero()
            || b*global::E.zero() != global::E.zero())
            err++;
    }
    end_test(err);
}

void Extension_test::test_even_tau()
{
    cout << "test even tau" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element sigma = global::E.random();
        Extension_element v = global::E.random();
        if (sigma.is_even() || v.is_even())
            /* we get here with probability (0.5)^(d-1) */
            continue;
        Extension_element e = v - sigma * util::tau(sigma, v);
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}

void Extension_test::test_is_even()
{
    cout << "test is even" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element e(0x0, global::randgen() & global::E.get_mask());
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}
