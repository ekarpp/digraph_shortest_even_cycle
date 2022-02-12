#include <iostream>

#include "../src/global.hh"
#include "gf_test.hh"

util::rand64bit global::randgen;

using namespace std;

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        cout << "pls only argument, n for GF(2^n)" << endl;
        return -1;
    }

    global::randgen.init(10);
    GF_test t(stoi(argv[1]));

    t.test_add_inverse();
    t.test_associativity();
    t.test_mul_inverse();
    t.test_mul_id();
    return 0;
}
