/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <vector>

#include "global.hh"
#include "solver.hh"
#include "graph.hh"
#include "gf.hh"
#include "util.hh"
#include "polynomial.hh"

using namespace std;

/* returns the length of the shortest even cycle in G.
 * if no even cycle exists, returns -1 */
int Solver::shortest_even_cycle(Graph G)
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
        delta[l] = G.get_A().pcc(gamma[l]);

    Polynomial p = util::poly_interpolation(gamma, delta);

    for (int k = 2; k <= G.get_n(); k += 2)
        if (p[G.get_n() - k] != global::F.zero())
            return k;

    return -1;
}
