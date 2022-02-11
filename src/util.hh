#ifndef UTIL_H
#define UTIL_H

#include <bitset>

using namespace std;

namespace util
{
    unsigned long long irred_poly(int deg);
    bool gcd1(int i, bitset<64> p);
}
#endif
