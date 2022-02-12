#include <iostream>
#include <vector>
#include <cmath>

#include "src/graph.hh"
#include "src/util.hh"
#include "src/global.hh"

using namespace std;

util::rand64bit global::randgen;

int main(int argc, char **argv)
{
    global::randgen.init(10);
    vector<vector<int>> vec = {
        {1},
        {0}
    };

    Graph G(vec);
    int d = 5 * ceil(log(G.get_n()) / log(2));
    if (d >= 64)
        return -1;

    unsigned long long p = util::irred_poly(d);
    return 0;
}
