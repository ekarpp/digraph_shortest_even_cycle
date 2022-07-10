/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "extension_test.hh"
#include "../../src/extension.hh"
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

        extension_repr ref = global::E.ref_mul(a.get_repr(), b.get_repr());
        extension_repr fast = global::E.fast_mul(a.get_repr(), b.get_repr());

        if (fast.hi != ref.hi || fast.lo != ref.lo)
            err++;
    }
    end_test(err);
}

void Extension_test::test_kronecker_mul()
{
    cout << "test kronecker mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        extension_repr ref = global::E.ref_mul(a.get_repr(), b.get_repr());
        extension_repr kron = global::E.kronecker_mul(a.get_repr(), b.get_repr());

        if (kron.hi != ref.hi || kron.lo != ref.lo)
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
        extension_repr v = global::E.fast_mul(a.get_repr(), b.get_repr());

        extension_repr euclid = global::E.euclid_rem(v);
        extension_repr intel = global::E.intel_rem(v);

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
        extension_repr a = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        extension_repr b = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        extension_repr mont = global::E.mont_reduce(
            global::E.mont_rem(
                global::E.mul(
                    global::E.mont_form(a),
                    global::E.mont_form(b)
                )
            )
        );

        extension_repr ref = global::E.euclid_rem(
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

#if GF2_bits == 16
void Extension_test::test_packed_intel_rem()
{
    cout << "test packed intel rem: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a;
        Extension_element b;
        extension_repr v = { 0, 0 };
        extension_repr euclid = { 0, 0 };

        for (int j = 0; j < 2; j++)
        {
            a = global::E.random();
            b = global::E.random();
            extension_repr tmp =
                global::E.fast_mul(a.get_repr(), b.get_repr());
            v.hi |= tmp.hi << (32*j);
            v.lo |= tmp.lo << (32*j);
            tmp = global::E.euclid_rem(tmp);
            euclid.hi |= tmp.hi << (32*j);
            euclid.lo |= tmp.lo << (32*j);
        }

        extension_repr intel = global::E.packed_intel_rem(v);

        if (euclid.hi != intel.hi || euclid.lo != intel.lo)
            err++;
    }
    end_test(err);
}

void Extension_test::test_packed_fast_mul()
{
    cout << "test packed fast mul: ";
    int err = 0;
    uint64_t mask = 0xFFFF;
    for (int i = 0; i < this->tests; i++)
    {
        extension_repr a;
        extension_repr b;
        extension_repr pa = { 0, 0 };
        extension_repr pb = { 0, 0 };
        extension_repr ref = { 0, 0 };

        for (int j = 0; j < 2; j++)
        {
            uint64_t r = global::randgen();

            a = {
                r & mask,
                (r >> 16) & mask
            };

            b = {
                (r >> 32) & mask,
                (r >> 48) & mask
            };

            pa.hi |= a.hi << (32*j);
            pa.lo |= a.lo << (32*j);

            pb.hi |= b.hi << (32*j);
            pb.lo |= b.lo << (32*j);

            extension_repr tmp = global::E.fast_mul(a, b);
            ref.hi |= tmp.hi << (32*j);
            ref.lo |= tmp.lo << (32*j);
        }

        extension_repr fast = global::E.packed_fast_mul(pa, pb);

        if (ref.hi != fast.hi || ref.lo != fast.lo)
            err++;
    }
    end_test(err);
}
#endif
