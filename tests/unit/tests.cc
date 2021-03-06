/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <getopt.h>

#include "../../src/global.hh"
#include "../../src/util.hh"
#include "../../src/gf.hh"
#include "../../src/extension.hh"

#include "gf_test.hh"
#include "extension_test.hh"
#include "fmatrix_test.hh"
#include "util_test.hh"
#include "solver_test.hh"
#include "ematrix_test.hh"
#include "geng_test.hh"

util::rand64bit global::randgen;
Extension global::E;
GF2n global::F;
bool global::output = false;

using namespace std;

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        cout << "-e for extension tests" << endl;
        cout << "-g for GF tests" << endl;
        cout << "-f for FMatrix tests" << endl;
        cout << "-x for EMatrix tests" << endl;
        cout << "-u for util tests" << endl;
        cout << "-s for solver tests" << endl;
        cout << "-d $int dimension of matrix" << endl;
        cout << "-n $int degree of modulo polynomial" << endl;
        cout << "-t $int how many times random tests are done" << endl;
        cout << "-c run geng test. \"geng -q $n | directg -q | listg -aq \" has to be piped to this." << endl;
        cout << "-r $int for seed to pseudo random generator" << endl;
        return 0;
    }

    bool et = false;
    bool gft = false;
    bool fmt = false;
    bool emt = false;
    bool ut = false;
    bool st = false;
    bool geng = false;
    int dim = 10;
    int tests = 10000;
    int opt;
    uint64_t seed = time(nullptr);

#if GF2_bits == 0
    int n = 10;
    while ((opt = getopt(argc, argv, "cxsuegfmn:d:t:r:")) != -1)
#else
    while ((opt = getopt(argc, argv, "cxsuegfmd:t:r:")) != -1)
#endif
    {
        switch (opt)
        {
#if GF2_bits == 0
        case 'n':
            n = stoi(optarg);
            break;
#endif
        case 'r':
            seed = stoll(optarg);
            break;
        case 'c':
            geng = true;
            break;
        case 'x':
            emt = true;
            break;
        case 'u':
            ut = true;
            break;
        case 'd':
            dim = stoi(optarg);
            break;
        case 's':
            st = true;
            break;
        case 'e':
            et = true;
            break;
        case 'g':
            gft = true;
            break;
        case 'f':
            fmt = true;
            break;
        case 't':
            tests = stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
#if GF2_bits == 0
    uint64_t mod = util::irred_poly(n);
    global::F.init(n, mod);
    global::E.init(n, mod);
#else
    global::F.init();
    global::E.init();
#endif

    if (et)
        Extension_test e(tests);
    if (gft)
        GF_test f;
    if (fmt)
        FMatrix_test fm(dim, tests);
    if (ut)
        Util_test u(dim, tests);
    if (st)
        Solver_test s(dim, tests);
    if (emt)
        EMatrix_test em(tests);
    if (geng)
        Geng_test g;

    return 0;
}
