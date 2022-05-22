/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "extension_test.hh"
#include "../../src/extension.hh"
#include "../../src/util.hh"
#include "../../src/global.hh"

using namespace std;


Extension_test::Extension_test(int tests)
{
    cout << "-----------------" << endl;
    cout << "TESTING EXTENSION" << endl;
    cout << "-----------------" << endl;
    this->tests = tests;
    this->run();
}

void Extension_test::test_add_inverse()
{
    cout << "add inverse: ";
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
    cout << "test associativity: ";
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
    cout << "test mul: ";
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

void Extension_test::test_fast_mul()
{
    cout << "test fast mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        uint64_2_t ref = global::E.ref_mul(a.get_repr(), b.get_repr());
        uint64_2_t fast = global::E.fast_mul(a.get_repr(), b.get_repr());

        if (fast.hi != ref.hi || fast.lo != ref.lo)
            err++;
    }
    end_test(err);
}

void Extension_test::test_intel_rem()
{
    cout << "test intel rem: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();
        uint64_2_t v = global::E.fast_mul(a.get_repr(), b.get_repr());

        uint64_2_t euclid = global::E.euclid_rem(v);
        uint64_2_t intel = global::E.intel_rem(v);

        if (euclid.hi != intel.hi || euclid.lo != intel.lo)
            err++;
    }
    end_test(err);
}

void Extension_test::test_mont_rem()
{
    cout << "test montgomery multiplication: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        uint64_2_t a = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        uint64_2_t b = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        uint64_2_t mont = global::E.mont_reduce(
            global::E.mont_rem(
                global::E.mul(
                    global::E.mont_form(a),
                    global::E.mont_form(b)
                )
            )
        );

        uint64_2_t ref = global::E.euclid_rem(
            global::E.mul(a, b)
        );

        if (ref.hi != mont.hi || ref.lo != mont.lo)
            err++;
    }
    end_test(err);
}

void Extension_test::test_even_tau()
{
    cout << "test even tau: ";
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
    cout << "test is even: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element e(0x0, global::randgen() & global::E.get_mask());
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}
