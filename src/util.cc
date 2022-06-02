/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <bitset>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>

#include "util.hh"
#include "extension.hh"
#include "gf.hh"
#include "global.hh"
#include "polynomial.hh"

using namespace std;

namespace util
{

    /* returns floor(log_2(a)) - 1 */
    int log2(uint64_t a)
    {
        /*
        int deg;
         why deg < 64 here?
        for (deg = 0; a >> deg && deg < 64; deg++);
        */
        return 63 - __builtin_clzl(a);
    }

    /* given an adjacency list for an undirected graph,
     * directs it such that edges are made one way
     * with direction chosen uniformly at random. */
    void direct_undirected(vector<vector<int>> &adj)
    {
        for (int u = 0; u < (int) adj.size(); u++)
        {
            vector<int> nbors = adj[u];
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
                vector<int>::iterator pos = find(
                    adj[del].begin(),
                    adj[del].end(),
                    keep
                );
                adj[del].erase(pos);
            }
        }
        return;
    }

    /* la grange interpolation with gamma and delta
     * note that we are in characteristic 2 and thus
     * - = +. done with the formula (3.3) here:
     * https://doi.org/10.1137/S0036144502417715 */
    Polynomial poly_interpolation(
        const vector<GF_element> &gamma,
        const vector<GF_element> &delta
    )
    {
        // assert(gamma.size() == delta.size())
        // assert(n > 2)
        int n = gamma.size();
        Polynomial interp(n - 1);

        /* weights*/
        vector<GF_element> w(n, global::F.one());
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
        vector<GF_element> P(n+1, global::F.zero());
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


    /* returns n distinct random elements from
     * global::F. (use LSFR?) */
    vector<GF_element> distinct_elements(int n)
    {
        vector<GF_element> vec(n);
        set<uint64_t> have;
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

    Extension_element tau(Extension_element sigma, Extension_element v)
    {
        return sigma.project().inv().lift() * v;
    }

    /* Ben-Or's irreducible polynomial generator.
     * returns irreducible polynomial of degree deg
     * in Z2[x] encoded as a bitstring
     */
    uint64_t irred_poly(int deg)
    {
        // assert(deg > 2)
        uint64_t mask = 1ll << (deg + 1);
        mask--;

        while (true)
        {
            bitset<64> p(global::randgen() & mask);
            p[deg] = true;

            int i;
            for (i = 1; i <= deg >> 1; i++)
            {
                if (!gcd1(i, p))
                    i = deg;
            }

            if (i < deg)
                return p.to_ullong();
        }
    }

    /* is gcd of x^(2^i) - x and p one (in Z2[x])*/
    bool gcd1(int i, bitset<64> p)
    {
        /* use set to represent polynomials for conveninece of
         * the methods and to save space as we have degree 2^i */
        /* better to use circular linked list if its created
         * for polynomial */
        set<int64_t> r;
        r.insert(1 << i);
        r.insert(1);

        set<int64_t> rn;
        for (int j = 0; j < 64; j++)
        {
            if (p[j])
                rn.insert(j);
        }

        if (*r.rbegin() < *rn.rbegin())
            r.swap(rn);

        /* standard Euclid's algo with Euclidean division */
        while (rn.size() != 0)
        {
            int64_t deg_r = *r.rbegin();
            int64_t deg_rn = *rn.rbegin();

            while (deg_r >= deg_rn)
            {
                set<int64_t>::reverse_iterator it = rn.rbegin();
                for ( ; it != rn.rend(); it++)
                {
                    int64_t id = *it + deg_r - deg_rn;
                    if (r.count(id) == 1)
                        r.erase(id);
                    else
                        r.insert(id);
                }
                if (r.empty())
                    deg_r = -1;
                else
                    deg_r = *r.rbegin();

            }

            rn.swap(r);
        }

        return r.size() == 1 && r.count(0) == 1;
    }
}
