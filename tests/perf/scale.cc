/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

#include "../../src/global.hh"
#include "../../src/graph.hh"
#include "../../src/solver.hh"
#include "../../src/util.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;
bool global::output = false;

int main(void)
{
    uint64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
    global::F.init();
    global::E.init();

    vector<vector<int>> graph = {
        {1},
        {0}
    };

    Solver s;
    for (int n = 3; n <= 64; n++)
    {
        graph.push_back({0});
        graph[n-2][0] = n-1;

        Graph G(graph);
        chrono::steady_clock::time_point start =
            chrono::steady_clock::now();
        int c = s.shortest_even_cycle(G);
        chrono::steady_clock::time_point end =
            chrono::steady_clock::now();
        cout << "c = " << c << " in time: " <<
            (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;
    }

    return 0;
}
