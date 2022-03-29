/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

#include "../../src/global.hh"
#include "../../src/gf.hh"
#include "../../src/util.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;

int main(int argc, char **argv)
{
    uint64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    uint64_t t = stoi(argv[argc - 1]);
    int d = 31;
    uint64_t mod = util::irred_poly(d);
    global::F.init(d, mod);

    GF_element a = global::F.random();
    GF_element b = global::F.one();
    chrono::steady_clock::time_point start =
        chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        b *= a;
    chrono::steady_clock::time_point end =
        chrono::steady_clock::now();
    b.print();
    cout << t << " multiplications in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;

    start = chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        a = a.inv();
    end = chrono::steady_clock::now();
    a.print();
    cout << t << " inversions in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;

    return 0;
}
