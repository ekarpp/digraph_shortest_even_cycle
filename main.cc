#include <iostream>
#include <vector>
#include <cmath>

#include "src/graph.hh"
#include "src/gf.hh"
#include "src/util.hh"

using namespace std;

int main(int argc, char **argv)
{
    vector<vector<int>> vec = {
        {1},
        {0}
    };

    Graph G(vec);
    int d = 5 * ceil(log(G.get_n()) / log(2));
    if (d >= 64)
        return -1;

    unsigned long long p = ben_or(d);

    return 0;
}
