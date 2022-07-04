/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_H
#define EXTENSION_H

#include <iostream>
#include <stdint.h>
#include <immintrin.h>

#include "gf.hh"
#include "util.hh"
#include "global.hh"
#include "bitvectors.hh"

/* forward declare */
class GF_element;
class Extension_element;

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
struct extension_repr
{
    uint64_t hi;
    uint64_t lo;
};

#if GF2_bits == 16
typedef uint128_t kronecker_form;
#else
struct kronecker_form
{
/* MSB */
    uint256_t big;
/* 32bit LSB */
    uint64_t small;
};
#endif

/* Extension of GF(2^n) to the ring E(4^n).
 * If GF(2^n) = Z2 / <g_2> for irreducible polynomial
 * g_2 of degree n, then if g_4 is g_2 but coefficients
 * projected to Z4 we have E(4^n) = Z4 / <g_4>. */
class Extension
{
private:
    int n;
    uint64_t mod;
    uint64_t mask;

    extension_repr mod_ast;
    extension_repr q_plus;

    extension_repr n_prime;
    extension_repr r_squared;

    /* euclidean division, only used once during initialization.
     * b has to be monic for this to work */
    extension_repr quo(extension_repr a, extension_repr b) const
    {
        extension_repr q = { 0, 0 };
        int dega = std::max(util::log2(a.lo), util::log2(a.hi));
        int degb = std::max(util::log2(b.lo), util::log2(b.hi));

        while (dega >= degb)
        {

            extension_repr s = {
                (a.hi & (1ll << dega)) >> degb,
                (a.lo & (1ll << dega)) >> degb
            };

            q = this->add(q, s);
            a = this->subtract(a, this->fast_mul(s, b));

            dega = std::max(util::log2(a.lo), util::log2(a.hi));
        }

        return q;
    }

public:
    Extension() {}

#if GF2_bits == 0
    void init(const int n, const uint64_t mod)
#else
    void init()
#endif
    {
#if GF2_bits == 16
        this->n = GF2_bits;

        /* x^16 + x^5 + x^3 + x^2 +  1 */
        this->mod = 0x1002D;
        this->mask = 0xFFFF;

        // N' = x^15 + x^14 + 3x^12 + x^7 + 3x^5 + 3x^4 + x^3 + x^2 + 3
        this->n_prime = { 0x1031, 0xD0BD };
        this->r_squared = { 0x018C, 0x0451 };
#elif GF2_bits == 32
        this->n = GF2_bits;

        /* x^32 + x^7 + x^3 + x^2 + 1 */
        this->mod = 0x10000008D;
        this->mask = 0xFFFFFFFF;

        // N' = x^30 + 3x^29 + x^28 + x^27 + x^26 + x^25 + x^23 + 3x^22 + x^21 + 2x^20 + 2x^19 + 3x^18 + x^15 + 2x^14 + 3x^13 + 2x^12 + x^10 + 3x^9 + 2x^8 + 2x^5 + 3x^4 + x^3 + x^2 + 3
        this->n_prime = { 0x205C7331, 0x7EE4A61D };
        this->r_squared = { 0x000006AC, 0x00004051 };
#else
        this->n = n;
        this->mod = mod;
        this->mask = (1ll << this->n) - 1;

        this->mod_ast = { 0, this->mod & this->mask };
        this->q_plus = this->quo({0, 1ull << (2*this->n)} , { 0, this->mod });

        this->r_squared = {
            0,
            1ull << (this->n * 2)
        };
        this->r_squared = this->rem(this->r_squared);

        // deg == 2*n
        extension_repr r = { 0x0, 1ull << (this->n*2) };
        int N = 1 << this->n;
        N -= 1;
        N *= 2;

        // deg <= n-1
        extension_repr r_rem = this->rem(r);
        extension_repr r_prime = r_rem;

        N--;
        long idx = 1ll << (util::log2(N) - 1);
        while (idx > 1)
        {
            r_prime = this->rem(this->mul(r_prime, r_prime));
            if (N & idx)
                r_prime = this->rem(this->mul(r_prime, r_rem));
            idx >>= 1;
        }

        // deg <= 3n - 1, overflow when 3n - 1 > 64 <=> n > 65 / 4 ~ 21.667
        this->n_prime = this->fast_mul(r, r_prime);
        this->n_prime = this->quo(this->n_prime, { 0x0, this->mod });
#endif

        if (global::output)
        {
            std::cout << "initialized E(4^" << this->n << ") with modulus: ";
            for (int i = this->n; i >= 0; i--)
            {
                if ((this->mod >> i) & 1)
                    std::cout << "1";
                else
                    std::cout << "0";
            }
            std::cout << std::endl;
        }
    }

    Extension_element zero() const;
    Extension_element one() const;
    Extension_element random() const;

    extension_repr rem(extension_repr a) const
    {
        return intel_rem(a);
    }

    extension_repr euclid_rem(extension_repr a) const
    {
        while (a.lo > this->mask || a.hi > this->mask)
        {
            int shift = std::max(util::log2(a.lo), util::log2(a.hi));
            shift -= this->n;
            /* mod has coefficients modulo 2, thus its negation
             * is just it applied to hi and lo (see negate function)*/
            a = this->add(a, { this->mod << shift, this->mod << shift });
        }
        return a;
    }

    /* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
    extension_repr intel_rem(extension_repr a) const
    {
        extension_repr hi = { a.hi >> this->n, a.lo >> this->n };
        extension_repr lo = { a.hi & this->mask, a.lo & this->mask };

#if GF2_bits == 16
        extension_repr tmp = { hi.hi >> 14, hi.lo >> 14 };
        tmp = this->add(tmp, { hi.hi >> 13, hi.lo >> 13 });
        tmp = this->add(tmp, { hi.hi >> 11, hi.lo >> 11 });
        tmp = this->subtract(hi, tmp);

        extension_repr r = this->add(tmp, { tmp.hi << 2, tmp.lo << 2 });
        r = this->add(r, { tmp.hi << 3, tmp.lo << 3 });
        r = this->add(r, { tmp.hi << 5, tmp.lo << 5 });
#elif GF2_bits == 32
        extension_repr tmp = { hi.hi >> 30, hi.lo >> 30 };
        tmp = this->add(tmp, { hi.hi >> 29, hi.lo >> 29 });
        tmp = this->add(tmp, { hi.hi >> 25, hi.lo >> 25 });
        tmp = this->subtract(hi, tmp);

        extension_repr r = this->add(tmp, { tmp.hi << 2, tmp.lo << 2 });
        r = this->add(r, { tmp.hi << 3, tmp.lo << 3 });
        r = this->add(r, { tmp.hi << 7, tmp.lo << 7 });
#else
        /* deg n-2 * deg n*/
        extension_repr r = this->mul(hi, this->q_plus);
        r = { r.hi >> this->n, r.lo >> this->n };
        /* deg n-1 * deg n - 2*/
        r = this->mul(r, this->mod_ast);
#endif
        r = { r.hi & this->mask, r.lo & this->mask };
        return this->subtract(lo, r);
    }

    extension_repr mont_rem(extension_repr a) const
    {
        /* d-1 deg * d-1 deg */
        extension_repr u = this->mul(a, this->n_prime);

        u.hi &= this->mask;
        u.lo &= this->mask;
        /* d deg * d-1 deg */
        extension_repr c = this->add(a, this->fast_mul(u, {0, this->mod}));

        c.hi >>= this->n;
        c.lo >>= this->n;
        return c;
    }

    extension_repr mont_form(extension_repr a) const
    {
        return this->mont_rem(this->mul(a, this->r_squared));
    }

    extension_repr mont_reduce(extension_repr a) const
    {
        return this->mont_rem(this->mul(a, { 0, 1 }));
    }

    extension_repr add(extension_repr a, extension_repr b) const
    {
        uint64_t carry = a.lo & b.lo;
        return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
    }

    extension_repr negate(extension_repr a) const
    {
        return {
            a.lo ^ a.hi,
            a.lo
        };
    }

    extension_repr subtract(extension_repr a, extension_repr b) const
    {
        return this->add(a, this->negate(b));
    }

    extension_repr mul(extension_repr a, extension_repr b) const
    {
        return kronecker_mul(a, b);
    }

    extension_repr ref_mul(extension_repr a, extension_repr b) const
    {
        extension_repr c = { 0, 0 };

        for (int i = 0; i <= global::E.get_n(); i++)
        {
            extension_repr tmp = { 0, 0 };
            if ((a.hi >> i) & 1)
                tmp.hi = this->mask;
            if ((a.lo >> i) & 1)
                tmp.lo = this->mask;

            /* 2 bit carryless multiplier */
            extension_repr aib = {
                (b.lo & tmp.hi) ^ (b.hi & tmp.lo),
                b.lo & tmp.lo
            };

            aib.lo <<= i;
            aib.hi <<= i;
            c = global::E.add(c, aib);
        }

        return c;
    }

    extension_repr fast_mul(extension_repr a, extension_repr b) const
    {
        /* clean this up */
        __m128i aa = _mm_set_epi64x(a.hi, a.lo);
        __m128i bb = _mm_set_epi64x(b.hi, b.lo);

        __m128i alobhi = _mm_clmulepi64_si128(aa, bb, 0x01);
        __m128i ahiblo = _mm_clmulepi64_si128(aa, bb, 0x10);

        uint64_t hi1 = _mm_extract_epi64(ahiblo, 0x0);
        uint64_t hi2 = _mm_extract_epi64(alobhi, 0x0);

        uint64_t hi = 0;
        uint64_t lo = 0;

        /* handle product of lo and lo */
        #pragma GCC unroll 32
#if GF2_bits == 16
        for (int i = 0; i <= GF2_bits; i++)
#else
        for (int i = 0; i <= 32; i++)
#endif
        {
            if ((b.lo >> i)&1)
            {
                hi ^= (a.lo << i) & lo;
                lo ^= (a.lo << i);
            }
        }

        return { hi1 ^ hi2 ^ hi, lo };
    }

    /* only works if deg <= 15 for a AND b */
    extension_repr kronecker_mul(extension_repr a, extension_repr b) const
    {
        extension_repr ret;
        /* we use different representation of polynomials than before here.
         * each bit string can be split to sets of 2 bits where each set
         * corresponds to a coefficient modulo 4. */
        kronecker_form aa = this->kronecker_substitution(a);
        kronecker_form bb = this->kronecker_substitution(b);
#if GF2_bits == 16
        uint256_t prod = bit::mul_128bit(aa, bb);

        /* first store the interesting bits to a uint64_t,
         * that is the first two bits of each 8 bit limb.
         * it fits, as we have deg <= 15+15 and each coefficient
         * uses two bits. */
        uint64_t extmask = 0x0303030303030303ull;
        uint64_t tmp = 0;
        for (int i = 0; i < 4; i++)
            tmp |= _pext_u64(prod.words[i], extmask) << (16*i);

        /* extract the usual hi/lo representation */
        uint64_t hiextmask = 0xAAAAAAAAAAAAAAAAull;
        uint64_t loextmask = 0x5555555555555555ull;
        ret.lo = _pext_u64(tmp, loextmask);
        ret.hi = _pext_u64(tmp, hiextmask);
#else
        uint512_t ahbh = bit::mul_256bit(aa.big, bb.big);
        uint512_t ahbl = bit::mul_256bit_64bit(aa.big, bb.small);
        uint512_t albh = bit::mul_256bit_64bit(bb.big, aa.small);
        uint64_t albl = aa.small * bb.small;

        uint576_t prod = bit::add_576bit(
            bit::widen_512bits(ahbh),
            bit::pad_words(ahbl, 4)
        );

        prod = bit::add_576bit(
            prod,
            bit::pad_words(albh, 4)
        );

        /* can't overflow */
        prod.words[8] += albl;

        /* extract */

        /* append zero bit to MSB of each word
         * to make each word contain exactly 7 coefficients:
         * 7*9 + 1 = 63 + 1 = 64 */

        for (int i = 8; i > 0; i--)
        {
            prod.words[i] <<= i;
            prod.words[i] |= prod.words[i-1] >> (64 - i);
        }

        uint64_t tmp[3];
        tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;

        uint64_t extmask = 0x00C06030180C0603ull;
        for (int i = 0; i < 4; i++)
            tmp[0] |= _pext_u64(prod.words[i], extmask) << (14*i);

        for (int i = 4; i < 8; i++)
            tmp[1] |= _pext_u64(prod.words[i], extmask) << (14*(i-4));

        tmp[2] = _pext_u64(prod.words[8], extmask);

        uint64_t hiextmask = 0xAAAAAAAAAAAAAAAAull;
        uint64_t loextmask = 0x5555555555555555ull;
        ret.hi = 0; ret.lo = 0;
        for (int i = 0; i < 3; i++)
        {
            ret.hi |= _pext_u64(tmp[i], hiextmask) << (28*i);
            ret.lo |= _pext_u64(tmp[i], loextmask) << (28*i);
        }
#endif
        return ret;
    }

    kronecker_form kronecker_substitution(extension_repr x) const
    {
        /* combine lo and hi to single uint64_t
         * where 2 bits represent single coefficient.
         * the "more traditional" bit representation for polynomials */
        uint64_t extmask = 0x5555555555555555ull;
        uint64_t comb = _pdep_u64(x.lo, extmask);
        comb |= _pdep_u64(x.hi, extmask << 1);

        kronecker_form kron;

#if GF2_bits == 16
        /* contains the "polynomial" after kronecker substitution.
         * for us it is sufficient that each coefficient has 8 bits,
         * (see details in thesis) thus we need 16*8 = 128 bits
         * for the polynomial after substitution. */

        extmask = 0x0303030303030303ull;
        kron.words[0] = _pdep_u64(comb & 0xFFFF, extmask);
        kron.words[1] = _pdep_u64((comb >> 16) & 0xFFFF, extmask);
#else
        /* each coefficients takes 9 bits.
         * we have <= 32 coefficients. */

        /* mask has 2x ones 7x zeros repeating */
        extmask = 0x00C06030180C0603ull;
        for (int i = 0; i < 4; i++)
            kron.big.words[i] = _pdep_u64((comb >> (i*14)) & 0x3FFF, extmask);

        /* remove the MSB zero from each word
         * to make the bitstring continuous */
        for (int i = 0; i < 3; i++)
        {
            kron.big.words[i] >>= i;
            kron.big.words[i] |= kron.big.words[i+1] << (63 - i);
        }
        kron.small = _pdep_u64((comb >> 56) & 0xFF, extmask);
        kron.big.words[3] >>= 3;
        kron.big.words[3] |= kron.small << 60;
        kron.small >>= 4;
#endif

        return kron;
    }

    int get_n() const { return this->n; }
    uint64_t get_mod() const { return this->mod; }
    uint64_t get_mask() const { return this->mask; }
};


class Extension_element
{
private:
    extension_repr repr;

public:
    Extension_element() { };

    Extension_element(const uint64_t lo, const uint64_t hi)
    {
        this->repr = { hi, lo };
    }

    Extension_element(const extension_repr repr)
    {
        this->repr = repr;
    }

    Extension_element(const Extension_element& e)
    {
        this->repr = { e.get_hi(), e.get_lo() };
    }

    Extension_element operator+(const Extension_element &other) const
    {
        return Extension_element(
            global::E.add(this->repr, other.get_repr())
        );
    }

    Extension_element &operator+=(const Extension_element &other)
    {
        this->repr = global::E.add(this->repr, other.get_repr());
        return *this;
    }

    Extension_element operator-(const Extension_element &other) const
    {
        return Extension_element(
            global::E.subtract(this->repr, other.get_repr())
        );
    }

    Extension_element &operator-=(const Extension_element &other)
    {
        this->repr = global::E.subtract(this->repr, other.get_repr());
        return *this;
    }

    Extension_element operator*(const Extension_element &other) const
    {
        extension_repr prod = global::E.mul(this->repr, other.get_repr());
        return Extension_element(global::E.rem(prod));
    }

    Extension_element &operator*=(const Extension_element &other)
    {
        extension_repr prod = global::E.mul(this->repr, other.get_repr());
        this->repr = global::E.rem(prod);
        return *this;
    }

    bool operator==(const Extension_element &other) const
    {
        return this->repr.lo == other.get_lo() && this->repr.hi == other.get_hi();
    }

    bool is_even() const
    {
        return this->repr.lo == 0x0;
    }

    /* used only on elements that are multiplied by two
     * thus we can just move the hi to low
     * maybe even just return gf element straight away
     * as this (probably?) gets anyways done after div2 */
    /* modify instead of returning new? */
    Extension_element div2() const
    {
        return Extension_element(this->repr.hi, 0x0);
    }

    GF_element project() const;

    uint64_t get_lo() const { return this->repr.lo; }
    uint64_t get_hi() const { return this->repr.hi; }
    extension_repr get_repr() const { return this->repr; }

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr.lo = other.get_lo();
        this->repr.hi = other.get_hi();
        return *this;
    }

    bool operator!=(const Extension_element &other) const
    {
        return !(*this == other);
    }

    void print() const
    {
        for (int i = 8; i >= 0; i--)
        {
            uint64_t v = (this->repr.hi >> i) & 1;
            v <<= 1;
            v |= (this->repr.lo >> i) & 1;
            std::cout << v;
        }
        std::cout << " ";
    }
};

namespace util
{
    inline Extension_element tau(Extension_element sigma, Extension_element v)
    {
        return sigma.project().inv().lift() * v;
    }
}

#endif
