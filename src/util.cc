/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <bitset>
#include <cmath>
#include <set>
#include <vector>

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
}
