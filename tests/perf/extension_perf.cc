/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/util.hh"

using namespace std;

util::rand64bit global::randgen;
Extension global::E;
GF2n global::F;
bool global::output = false;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "-s $int for seed" << endl;
        cout << "-t $int for amount of tests" << endl;
        return 0;
    }

    uint64_t seed = time(nullptr);

    uint64_t t = 1;
    int opt;
    while ((opt = getopt(argc, argv, "s:t:")) != -1)
    {
        switch (opt)
        {
        case 's':
            seed = stoi(optarg);
            break;
        case 't':
            t = stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
#if GF2_bits == 0
    int n = 20;
    uint64_t mod = util::irred_poly(n);
    global::E.init(n, mod);
#else
    global::E.init();
#endif

    vector<uint64_2_t> a(t);
    vector<uint64_2_t> b(t);
    vector<uint64_2_t> p(t);
    vector<uint64_2_t> r(t);

    for (uint64_t i = 0; i < t; i++)
    {
        a[i] = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        b[i] = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };
    }

    double start;
    double end;
    double delta;
    double mhz;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.ref_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;
#if GF2_bits == 16
    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.kronecker_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;
#endif
    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.fast_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::E.mont_rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::E.euclid_rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::E.intel_rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " intel remainders in time " <<
        delta << " s or " << mhz << " Mhz" << endl;
    return 0;
}
