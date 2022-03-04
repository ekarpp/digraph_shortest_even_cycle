/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <getopt.h>

#include "../../src/global.hh"
#include "../../src/util.hh"
#include "gf_test.hh"
#include "extension_test.hh"
#include "matrix_test.hh"
#include "fmatrix_test.hh"
#include "util_test.hh"
#include "solver_test.hh"

util::rand64bit global::randgen;
Extension global::E;
GF2n global::F;

using namespace std;

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        cout << "-e for extension tests" << endl;
        cout << "-g for GF tests" << endl;
        cout << "-f for FMatrix tests" << endl;
        cout << "-m for Matrix tests" << endl;
        cout << "-u for util tests" << endl;
        cout << "-s for solver tests" << endl;
        cout << "-d $int dimension of matrix" << endl;
        cout << "-n $int degree of modulo polynomial" << endl;
        cout << "-t $int how many times random tests are done" << endl;
        return 0;
    }

    bool et = false;
    bool gft = false;
    bool mt = false;
    bool fmt = false;
    bool ut = false;
    bool st = false;
    int n = 10;
    int dim = 10;
    int tests = 10000;
    int opt;
    while ((opt = getopt(argc, argv, "suegfmd:n:t:")) != -1)
    {
        switch (opt)
        {
        case 'u':
            ut = true;
            break;
        case 'n':
            n = stoi(optarg);
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
        case 'm':
            mt = true;
            break;
        case 't':
            tests = stoi(optarg);
            break;
        }
    }

    uint64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    int64_t p = util::irred_poly(n);

    global::E.init(n, p);
    global::F.init(n, p);

    if (et)
        Extension_test e(n);
    if (gft)
        GF_test f(n);
    if (mt)
        Matrix_test m(dim, tests);
    if (fmt)
        FMatrix_test fm(dim, tests);
    if (ut)
        Util_test u(n, tests);
    if (st)
        Solver_test s(n, tests);

    return 0;
}
