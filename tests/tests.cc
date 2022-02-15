#include <iostream>

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
    if (argc != 2)
    {
        cout << "pls only argument, n for GF(2^n)" << endl;
        return -1;
    }
    int64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    int n = stoi(argv[1]);
    int64_t p = util::irred_poly(n);

    global::E.init(n, p);
    global::F.init(n, p);

    Extension_test e(n);
    GF_test f(n);
    Matrix_test m(5);
    FMatrix_test fm(5);

    return 0;
}
