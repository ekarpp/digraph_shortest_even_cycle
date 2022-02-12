#include <bitset>
#include <cmath>
#include <set>

#include "util.hh"

using namespace std;


namespace util
{
    /* Ben-Or's irreducible polynomial generator.
     * returns irreducible polynomial of degree deg
     * in Z2[x] encoded as a bitstring
     */
    int64_t irred_poly(int deg)
    {
        // assert(deg > 2)
        bitset<64> p;
        p[deg] = true;
        srand(1);

        while (true)
        {
            int i;
            for (i = 0; i < deg; i++)
                p[i] = rand() & 1;
            for (i = 1; i <= deg >> 1; i++)
            {
                if (!gcd1(i, p))
                    i = deg;
            }

            if (i < deg)
                break;
        }

        return p.to_ullong();
    }

    /* is gcd of x^(2^i) - x and p one (in Z2[x])*/
    bool gcd1(int i, bitset<64> p)
    {
        /* use set to represent polynomials for conveninece of
         * the methods and to save space as we have degree 2^i */
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

    /* returns r s.t. for some q,
     * a = q*b + r is the division equation (in Z(2^n))
     * bdeg is degree of b (can be computed, but is stored in field anyways)
     */
    // a needs to be 128 bit for support up to GF(2^64)
    // now just GF(2^32)
    int64_t modz2(int64_t a, int64_t b, int bdeg)
    {
        // assert(b != 0)
        int64_t mask = (1 << bdeg) - 1;
        while (a > mask)
        {
            int shift;
            for (shift = 0; a >> shift; shift++);
            shift -= 1 + bdeg;
            /* shift = deg(a) - deg(b) */
            a ^= (b << shift);
        }
        return a;
    }
}
