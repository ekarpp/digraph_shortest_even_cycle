/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <bitset>
#include <vector>

#include "extension.hh"
#include "gf.hh"
#include "polynomial.hh"

using namespace std;

namespace util
{
    int log2(uint64_t a);

    uint64_t irred_poly(int deg);
    bool gcd1(int i, bitset<64> p);

    std::vector<GF_element> distinct_elements(int n);

    Polynomial poly_interpolation(
        const vector<GF_element> &gamma,
        const vector<GF_element> &delta
    );

    Extension_element tau(Extension_element sigma, Extension_element v);

    class rand64bit
    {
    private:
        mt19937_64 gen;
        uniform_int_distribution<uint64_t> dist;
    public:
        rand64bit() {}
        void init(uint64_t seed) { gen = mt19937_64(seed); }
        uint64_t operator() () { return dist(gen); }
    };

}
#endif
