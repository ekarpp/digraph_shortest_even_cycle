/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <vector>
#include <iostream>

#include "global.hh"
#include "solver.hh"
#include "graph.hh"
#include "gf.hh"
#include "util.hh"
#include "polynomial.hh"

using namespace std;

/* returns the length of the shortest even cycle in G.
 * if no even cycle exists, returns -1 */
int Solver::shortest_even_cycle(Graph G, bool out)
{
    vector<GF_element> gamma(G.get_n() + 1);
    vector<GF_element> delta(G.get_n() + 1);

    uint64_t v = 0x1;

    for (int i = 0; i <= G.get_n(); i++)
    {
        gamma[i] = GF_element(v);
        /* lazy, fix later */
        v = global::F.rem(v + 3);
    }

    for (int l = 0; l <= G.get_n(); l++)
    {
        delta[l] = G.get_A().pcc(gamma[l]);
        if (out)
            cout << l+1 << "/" << G.get_n()+1 << endl;
    }

    Polynomial p = util::poly_interpolation(gamma, delta);

    for (int k = 2; k <= G.get_n(); k += 2)
        if (p[G.get_n() - k] != global::F.zero())
            return k;

    return -1;
}

/* brute force approach with DFS to find the
 * shortest even cycle in G. returns the length
 * of the shortest even cycle or -1 if no such
 * cycle exists */
int Solver::shortest_even_cycle_brute(Graph G)
{
    int len = G.get_n() + 1;

    for (int v = 0; v < G.get_n(); v++)
    {
        vector<bool> visited(G.get_n(), false);
        G.dfs_cycle(v, 1, v, visited, &len);
    }

    return (len == G.get_n() + 1) ? -1 : len;
}
