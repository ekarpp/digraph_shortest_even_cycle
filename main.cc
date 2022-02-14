#include <iostream>
#include <vector>
#include <cmath>

#include "src/gf.hh"
#include "src/global.hh"
#include "src/graph.hh"
#include "src/util.hh"
#include "src/fmatrix.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;

int main(int argc, char **argv)
{
    global::randgen.init(10);
    vector<vector<int>> vec = {
        {1},
        {0}
    };

    Graph G(vec);
    int d = 5 * ceil(log(G.get_n()) / log(2));
    /* modz2 and gfmul only support upto 32 */
    if (d > 32)
        return -1;

    int64_t poly = util::irred_poly(d);
    global::F.init(d, poly);
    global::E.init(d, poly);

    return 0;
}
