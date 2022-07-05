/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/gf.hh"

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
        cout << "-t $int for amount of tests, exponent for two" << endl;
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
            t <<= stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
#if GF2_bits == 0
    int n = 24;
    cout << "E^" << n << endl;
    uint64_t mod = util::irred_poly(n);
    global::E.init(n, mod);
#else
    global::E.init();
#endif

    vector<extension_repr> a(t);
    vector<extension_repr> b(t);
    vector<extension_repr> p(t);
    vector<extension_repr> r(t);

    double start;
    double end;

    start = omp_get_wtime();
#if GF2_bits == 16
    #define REPEATS 2
#else
    #define REPEATS 1
#endif
    #pragma omp parallel for
    for (uint64_t i = 0; i < t; i += REPEATS)
    {
        uint64_t aa = global::randgen();
        uint64_t bb = global::randgen();
        for (uint64_t j = 0; j < REPEATS; j++)
        {
            uint64_t alo = aa >> (2*GF2_bits*j);
            uint64_t ahi = aa >> (GF2_bits*(2*j+1));
            a[i+j] = {
                ahi & global::E.get_mask(),
                alo & global::E.get_mask()
            };

            uint64_t blo = bb >> (2*GF2_bits*j);
            uint64_t bhi = bb >> (GF2_bits*(2*j+1));
            b[i+j] = {
                bhi & global::E.get_mask(),
                blo & global::E.get_mask()
            };
        }
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    double delta;
    double mhz;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.ref_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.kronecker_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::E.fast_mul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::E.mont_rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::E.euclid_rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
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
