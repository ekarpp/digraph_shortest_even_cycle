/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>
#include <immintrin.h>

#include "../../src/global.hh"
#include "../../src/gf.hh"
#include "../../src/extension.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;
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
    int n = 10;
    uint64_t mod = util::irred_poly(n);
    global::F.init(n, mod);
#else
    global::F.init();
#endif

    double start, end, delta, mhz;

#if GF2_bits == 16
#define WIDTH 16
    vector<__m512i> aw(t);
    vector<__m512i> bw(t);

    for (uint64_t i = 0; i < t; i++)
    {
        uint64_t va[WIDTH/2];
        uint64_t vb[WIDTH/2];

        uint64_t packmask = (0xFFFFull << 32) | 0xFFFF;
        for (uint64_t j = 0; j < WIDTH / 2; j++)
        {
            va[j] = global::randgen() & packmask;
            vb[j] = global::randgen() & packmask;
        }

        aw[i] = _mm512_set_epi64(va[0], va[1], va[2], va[3], va[4], va[5], va[6], va[7]);
        bw[i] = _mm512_set_epi64(vb[0], vb[1], vb[2], vb[3], vb[4], vb[5], vb[6], vb[7]);
    }


    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        aw[i] = global::F.wide_mul(aw[i], bw[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = (WIDTH*t) / delta;
    mhz /= 1e6;

    cout << t*WIDTH << " wide multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

#endif

    vector<uint64_t> a(t);
    vector<uint64_t> b(t);
    vector<uint64_t> p(t);
    vector<uint64_t> r(t);

    vector<GF_element> aa(t);
    vector<GF_element> bb(t);

    for (uint64_t i = 0; i < t; i++)
    {
        a[i] = global::randgen() & global::F.get_mask();
        b[i] = global::randgen() & global::F.get_mask();
        aa[i] = global::F.random();
        bb[i] = global::F.random();
    }

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        aa[i] *= bb[i];
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " multiplications (whole) in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::F.clmul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F.rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " remainder in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
#ifdef PARALLEL
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F.ext_euclid(r[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " inversion in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    return 0;
}
