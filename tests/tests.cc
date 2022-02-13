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
    e.run();

    GF_test f(stoi(argv[1]));
    f.run();

    return 0;
}
