#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <bitset>

#include "extension.hh"

using namespace std;

namespace util
{
    uint64_t irred_poly(int deg);
    bool gcd1(int i, bitset<64> p);

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
