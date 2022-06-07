/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <set>
#include <bitset>
#include <vector>

#include "extension.hh"
#include "gf.hh"
#include "polynomial.hh"

namespace util
{
/*    inline int log2(uint64_t a)
    {
        return 63 - __builtin_clzl(a);
    }*/

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

    /* returns n distinct random elements from
     * global::F. (use LSFR?) */
    inline std::vector<GF_element> distinct_elements(int n)
    {
        std::vector<GF_element> vec(n);
        std::set<uint64_t> have;
        for (int i = 0; i < n; i++)
        {
            GF_element e = global::F.random();
            while (have.count(e.get_repr()) == 1)
                e = global::F.random();
            vec[i] = e;
            have.insert(e.get_repr());
        }
        return vec;
    }

    /* la grange interpolation with gamma and delta
     * note that we are in characteristic 2 and thus
     * - = +. done with the formula (3.3) here:
     * https://doi.org/10.1137/S0036144502417715 */
    inline Polynomial poly_interpolation(
        const std::vector<GF_element> &gamma,
        const std::vector<GF_element> &delta
    )
    {
        // assert(gamma.size() == delta.size())
        // assert(n > 2)
        int n = gamma.size();
        Polynomial interp(n - 1);

        /* weights*/
        std::vector<GF_element> w(n, global::F.one());
        for (int j = 1; j < n; j++)
        {
            for (int k = 0; k < j; k++)
            {
                w[k] *= gamma[k] + gamma[j];
                w[j] *= gamma[k] + gamma[j];
            }
        }

        for (int j = 0; j < n; j++)
            w[j].inv_in_place();

        /* main polynomial [ prod_{i} (x + gamma_i) ]*/
        std::vector<GF_element> P(n+1, global::F.zero());
        P[n] = global::F.one();
        P[n-1] = gamma[0];
        for (int i = 1; i < n; i++)
        {
            for (int j = n - i - 1; j < n - 1; j++)
                P[j] += gamma[i] * P[j+1];
            P[n - 1] += gamma[i];
        }

        for (int i = 0; i < n; i++)
        {
            Polynomial tmp(P);
            tmp.div(gamma[i]);
            tmp *= w[i] * delta[i];
            interp += tmp;
        }

        return interp;
    }

    inline Extension_element tau(Extension_element sigma, Extension_element v)
    {
        return sigma.project().inv().lift() * v;
    }

    uint64_t irred_poly(int deg);
    bool gcd1(int i, std::bitset<64> p);
}
#endif
