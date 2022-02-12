#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <bitset>

using namespace std;

namespace util
{
    int64_t irred_poly(int deg);
    bool gcd1(int i, bitset<64> p);


    class rand64bit
    {
    private:
        mt19937_64 gen;
        uniform_int_distribution<uint64_t> dist;
    public:
        rand64bit() {}
        void init(int seed) { gen = mt19937_64(seed); }
        int64_t operator() () { return dist(gen); }
    };

}
#endif
