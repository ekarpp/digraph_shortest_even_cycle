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
    /* la grange interpolation with gamma and delta
     * note that we are in characteristic 2 and thus
     * - = + */
    Polynomial poly_interpolation(
        const vector<GF_element> &gamma,
        const vector<GF_element> &delta
    )
    {
        // assert(gamma.size() == delta.size())
        // assert(n > 2)
        int n = gamma.size();
        Polynomial ret(n - 1);

        for (int i = 0; i < n; i++)
        {

            int start = (i == 0) ? 1 : 0;
            int deg = n - 1;
            /* init poly with coeff 1 deg n-1 */
            Polynomial tmp(deg);
            tmp(deg, global::F.one());
            deg--;

            GF_element coeff = gamma[start];
            GF_element prod = gamma[start];
            GF_element quotient = gamma[i] + gamma[start];
            start++;
            if (i == start)
                start++;
            coeff += gamma[start];
            tmp(deg, coeff);
            deg--;

            prod *= gamma[start];
            quotient *= gamma[i] + gamma[start];

            for (int j = start + 1; j < n; j++)
            {
                if (i == j)
                    continue;

                coeff *= gamma[j];
                tmp(deg, coeff);
                deg--;
                quotient *= gamma[i] + gamma[j];
                prod *= gamma[j];
            }
            /* add to poly prod with deg 0 */
            tmp(0, prod);

            tmp *= delta[i] / quotient;
            ret += tmp;
        }

        return ret;
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

    Extension_element tau(Extension_element sigma, Extension_element v)
    {
        return sigma.project().inv().lift() * v;
    }
}
