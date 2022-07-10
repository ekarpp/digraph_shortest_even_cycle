/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/gf.hh"

#ifndef PARALLEL
#define PARALLEL 0
#endif

#if GF2_bits == 16
#define VAL_PER_WORD 4
#else
#define VAL_PER_WORD 2
#endif

#if GF2_bits == 0
#define SHIFT_N (n)
#else
#define SHIFT_N GF2_bits
#endif

#define WARMUP (1 << 20)

#define EXTRACT_MUL_LOOP(op)                            \
{                                                       \
    for (int j = 0; j < VAL_PER_WORD; j++)              \
    {                                                   \
        extension_repr lhs =                            \
            a[i].shiftr_and(                            \
                j*SHIFT_N,                              \
                global::E.get_mask()                    \
            );                                          \
        extension_repr rhs =                            \
            b[i].shiftr_and(                            \
                j*SHIFT_N,                              \
                global::E.get_mask()                    \
            );                                          \
        op(rhs, lhs);                                   \
    }                                                   \
}

#define EXTRACT_REM_LOOP(op)                            \
{                                                       \
    for (int j = 0; j < VAL_PER_WORD/2; j++)            \
    {                                                   \
        extension_repr lhs =                            \
            a[i].shiftr_and(                            \
                j*SHIFT_N*2,                            \
                global::E.get_mask() |                  \
                (global::E.get_mask() << SHIFT_N)       \
            );                                          \
        extension_repr rhs =                            \
            b[i].shiftr_and(                            \
                j*SHIFT_N*2,                            \
                global::E.get_mask() |                  \
                (global::E.get_mask() << SHIFT_N)       \
            );                                          \
        op(rhs);                                        \
        op(lhs);                                        \
    }                                                   \
}

#define BENCH_MUL(mul_func)                                 \
{                                                           \
    extension_repr w = {0, 0};                              \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            EXTRACT_MUL_LOOP(w = mul_func)                  \
    } else {                                                \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            EXTRACT_MUL_LOOP(w = mul_func)                  \
    }                                                       \
    start = omp_get_wtime();                                \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < t; i++)                    \
            EXTRACT_MUL_LOOP(aa[i+j] = mul_func)            \
    } else {                                                \
        for (uint64_t i = 0; i < t; i++)                    \
            EXTRACT_MUL_LOOP(aa[i+j] = mul_func)            \
    }                                                       \
    end = omp_get_wtime();                                  \
    if (exp == 0x6fabc73829101) {                           \
        cout << aa[exp].hi << aa[exp].lo << endl;           \
        cout << w.hi << w.lo << endl;                       \
    }                                                       \
    delta = (end - start);                                  \
    mhz = t / delta;                                        \
    mhz /= 1e6;                                             \
}

#define BENCH_REM(rem_func)                                 \
{                                                           \
    extension_repr w = {0, 0};                              \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            EXTRACT_REM_LOOP(w = rem_func)                  \
    } else {                                                \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            EXTRACT_REM_LOOP(w = rem_func)                  \
    }                                                       \
    start = omp_get_wtime();                                \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < t; i++)                    \
            EXTRACT_REM_LOOP(aa[i+j] = rem_func)            \
    } else {                                                \
        for (uint64_t i = 0; i < t; i++)                    \
            EXTRACT_REM_LOOP(aa[i+j] = rem_func)            \
    }                                                       \
    end = omp_get_wtime();                                  \
    if (exp == 0x6fabc73829101) {                           \
        cout << aa[exp].hi << aa[exp].lo << endl;           \
        cout << w.hi << w.lo << endl;                       \
    }                                                       \
    delta = (end - start);                                  \
    mhz = t / delta;                                        \
    mhz /= 1e6;                                             \
}

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
    int exp = 1;
    int opt;
    while ((opt = getopt(argc, argv, "s:t:")) != -1)
    {
        switch (opt)
        {
        case 's':
            seed = stoi(optarg);
            break;
        case 't':
            exp = stoi(optarg);
            t <<= exp;
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

    t /= VAL_PER_WORD;

    vector<extension_repr> a(t);
    vector<extension_repr> b(t);
    vector<extension_repr> aa(t*VAL_PER_WORD);

    double start;
    double end;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
    {
        a[i] = {
            global::randgen(),
            global::randgen()
        };
        aa[i] = a[i];
        b[i] = {
            global::randgen(),
            global::randgen()
        };
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    double delta;
    double mhz;

    BENCH_MUL(global::E.ref_mul);

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E.kronecker_mul);

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E.fast_mul);

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E.mont_rem);

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E.euclid_rem);

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E.intel_rem);

    cout << t << " intel remainders in time " <<
        delta << " s or " << mhz << " Mhz" << endl;
    return 0;
#if GF2_bits == 16
    uint64_t pack_mask = 0xFFFF | (0xFFFFull << 32);
    t /= 2;
    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
    {
        uint64_t ar = global::randgen();
        uint64_t br = global::randgen();

        a[i] = {
            ar & pack_mask,
            (ar >> 16) & pack_mask
        };
        aa[i] = a[i];
        b[i] = {
            br & pack_mask,
            (br >> 16) & pack_mask
        };
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    BENCH_MUL(global::E.packed_fast_mul);

    cout << 2*t << " pack fast multiplications in time: " <<
        delta << " s or " << 2*mhz << " Mhz" << endl;

    BENCH_REM(global::E.packed_intel_rem);

    cout << 2*t << " pack intel remainders in time " <<
        delta << " s or " << 2*mhz << " Mhz" << endl;
#endif

    return 0;
}
