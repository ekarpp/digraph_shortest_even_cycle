/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "gf_test.hh"
#include "../../src/gf.hh"
#include "../../src/global.hh"
#include "../../src/extension.hh"

using namespace std;

GF_test::GF_test()
{
    cout << "----------" << endl;
    cout << "TESTING GF" << endl;
    cout << "----------" << endl;

    this->run();
}

void GF_test::test_add_inverse()
{
    cout << "add inverse: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e + e != global::F.zero()
            || e - e != global::F.zero())
            err++;
        i++;
    }
    end_test(err);
}

void GF_test::test_associativity()
{
    cout << "test associativity: ";
    int err = 0;
    for (int i = 0; i < 10000; i++)
    {
        GF_element a = global::F.random();
        GF_element b = global::F.random();
        GF_element c = global::F.random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    this->end_test(err);
}

void GF_test::test_mul_id()
{
    cout << "mul with id: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e * global::F.one() != e)
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_mul_inverse()
{
    cout << "mul with inverse: ";
    int err = 0;
    uint64_t i = 1;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e / e != global::F.one())
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_lift_project()
{
    cout << "lift project: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        Extension_element b(i, global::randgen() & global::E.get_mask());
        Extension_element c(i, 0b0);
        if (e.lift().project() != e
            || b.project() != e
            || e.lift() != c)
            err++;
        i++;
    }
    this->end_test(err);
}
