/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <bitset>
#include <vector>

#include "global.hh"

namespace util
{
    inline int log2(uint64_t a)
    {
        return 63 - __builtin_clzl(a);
    }

    /* given an adjacency list for an undirected graph,
     * directs it such that edges are made one way
     * with direction chosen uniformly at random. */
    inline void direct_undirected(std::vector<std::vector<int>> &adj)
    {
        for (int u = 0; u < (int) adj.size(); u++)
        {
            std::vector<int> nbors = adj[u];
            for (int i = 0; i < (int) nbors.size(); i++)
            {
                int v = nbors[i];
                if (u > v)
                    continue;
                /* rndom choose del */
                int keep = u;
                int del = v;
                if (global::randgen() % 2)
                {
                    keep = v;
                    del = u;
                }
                std::vector<int>::iterator pos = std::find(
                    adj[del].begin(),
                    adj[del].end(),
                    keep
                );
                adj[del].erase(pos);
            }
        }
        return;
    }

    uint64_t irred_poly(int deg);
    bool gcd1(int i, std::bitset<64> p);
}
#endif
