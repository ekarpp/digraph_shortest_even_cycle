/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <getopt.h>

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
    while ((opt = getopt(argc, argv, "mifs:t:")) != -1)
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
    global::E.init();

    vector<uint64_2_t> a(t);
    vector<uint64_2_t> b(t);
    vector<uint64_2_t> c(t);

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

    chrono::steady_clock::time_point start;
    chrono::steady_clock::time_point end;

    start = chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        c[i] = global::E.euclid_rem(global::E.ref_mul(a[i], b[i]));
    end = chrono::steady_clock::now();
    cout << t << " multiplications (euclid + ref) in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;


    start = chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        c[i] = global::E.intel_rem(global::E.fast_mul(a[i], b[i]));
    end = chrono::steady_clock::now();

    cout << t << " multiplications (intel + fast) in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;


    start = chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        c[i] = global::E.mont_mul(a[i], b[i]);
    end = chrono::steady_clock::now();

    cout << t << " multiplications (montgomery + fast) in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;
    return 0;
}
