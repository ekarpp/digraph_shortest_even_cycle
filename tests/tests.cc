#include <iostream>

#include "../src/global.hh"
#include "../src/util.hh"
#include "gf_test.hh"
#include "extension_test.hh"

util::rand64bit global::randgen;

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

    Extension_test e(stoi(argv[1]));
    e.test_add_inverse();
    e.test_associativity();
    e.test_mul();

    GF_test f(stoi(argv[1]));
    f.test_add_inverse();
    f.test_associativity();
    f.test_mul_inverse();
    f.test_mul_id();
    return 0;
}
