/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>

#include "util_test.hh"
#include "../src/gf.hh"
#include "../src/polynomial.hh"
#include "../src/util.hh"
#include "../src/global.hh"

using namespace std;

Util_test::Util_test(int tests = 0)
{
    if (tests)
        this->tests = tests;

    cout << "------------" << endl;
    cout << "TESTING UTIL" << endl;
    cout << "------------" << endl;

    this->run();
}

void Util_test::test_interpolation()
{
    cout << "polynomial interpolation: ";
    int err = 0;
    int n = 3;
    for (int t = 0; t < this->tests; t++)
    {
        uint64_t g = global::F.rem(global::randgen());
        uint64_t d = global::F.rem(global::randgen());

        vector<GF_element> gamma(n);
        vector<GF_element> delta(n);

        for (int i = 0; i < n; i++)
        {
            gamma[i] = GF_element(g);
            g = global::F.rem(g + 3);

            delta[i] = GF_element(d);
            d = global::F.rem(d + 3);
        }

        Polynomial p = util::poly_interpolation(gamma, delta);

        for (int i = 0; i < n; i++)
        {
            if (p.eval(gamma[i]) != delta[i])
            {
                err++;
                break;
            }
        }
    }
    end_test(err);
}
