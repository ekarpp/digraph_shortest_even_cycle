#include <iostream>
#include <getopt.h>

#include "../src/global.hh"
#include "../src/util.hh"
#include "gf_test.hh"
#include "extension_test.hh"
#include "matrix_test.hh"
#include "fmatrix_test.hh"

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
        cout << "-d $int dimension of matrix" << endl;
        cout << "-n $int degree of modulo polynomial" << endl;
        cout << "-t $int how many times random tests are done" << endl;
        return 0;
    }

    bool et = false;
    bool gft = false;
    bool mt = false;
    bool fmt = false;
    int n = 10;
    int dim = 10;
    int tests = 10000;
    int opt;
    while ((opt = getopt(argc, argv, "egfmd:n:t:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            n = stoi(optarg);
            break;
        case 'd':
            dim = stoi(optarg);
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

    return 0;
}
