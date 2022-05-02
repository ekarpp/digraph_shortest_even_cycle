/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <bitset>
#include <vector>

#include "extension.hh"
#include "gf.hh"
#include "polynomial.hh"

namespace util
{
    int log2(uint64_t a);

    std::vector<GF_element> distinct_elements(int n);

    Polynomial poly_interpolation(
        const std::vector<GF_element> &gamma,
        const std::vector<GF_element> &delta
    );

    Extension_element tau(Extension_element sigma, Extension_element v);

    uint64_t irred_poly(int deg);
    bool gcd1(int i, std::bitset<64> p);
}
#endif
