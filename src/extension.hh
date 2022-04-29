/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_H
#define EXTENSION_H

#include <iostream>
#include <stdint.h>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"

/* forward declare */
class GF_element;
class Extension_element;

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
/* vectorize with AVX? NO */
struct uint64_2_t
{
    uint64_t hi;
    uint64_t lo;
};

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

    uint64_2_t mod_ast;
    uint64_2_t q_plus;

    uint64_2_t n_prime;

    /* euclidean division, only used once during initialization.
     * b has to be monic for this to work */
    uint64_2_t quo(uint64_2_t a, uint64_2_t b) const
    {
        uint64_2_t q = { 0, 0 };
        int degb = this->n;
        int dega = 63 - std::min(__builtin_clzl(a.lo), __builtin_clzl(a.hi));
        while (dega >= degb)
        {
            uint64_2_t s = {
                (a.hi & (1ll << dega)) >> degb,
                (a.lo & (1ll << dega)) >> degb
            };

            q = this->add(q, s);
            a = this->add(a, this->negate(this->mul(b, s)));

            dega = 63 - std::min(__builtin_clzl(a.lo), __builtin_clzl(a.hi));
        }

        return q;
    }

public:
    Extension() {}

    void init()
    {
        this->n = GF2_bits;
#if GF2_bits == 16
        /* x^16 + x^5 + x^3 + x^2 +  1 */
        this->mod = 0x1002D;
        this->mask = 0xFFFF;
        this->q_plus = { 0, 0x1002D };
        this->mod_ast = { 0, 0x2D };
        // N' = x^15 + x^14 + 3x^12 + x^7 + 3x^5 + 3x^4 + x^3 + x^2 + 3
        this->n_prime = { 0x1031, 0xD0BD };
#elif GF2_bits == 32
        /* x^32 + x^7 + x^3 + x^2 + 1 */
        this->mod = 0x10000008D;
        this->mask = 0xFFFFFFFF;
        this->q_plus = { 0, 0x10000008D };
        this->mod_ast = { 0, 0x8D };
        this->n_prime = { 0x1031, 0xD0BD };
#else
        GF2_bits_eq_16_or_32
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

    uint64_2_t rem(uint64_2_t a) const
    {
        return intel_rem(a);
    }

    uint64_2_t euclid_rem(uint64_2_t a) const
    {
        while (a.lo > this->mask || a.hi > this->mask)
        {
            int shift = 63 - std::min(__builtin_clzl(a.lo), __builtin_clzl(a.hi));
            shift -= this->n;
            /* mod has coefficients modulo 2, thus its negation
             * is just it applied to hi and lo (see negate function)*/
            a = this->add(a, { this->mod << shift, this->mod << shift });
        }
        return a;
    }

    /* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
    uint64_2_t intel_rem(uint64_2_t a) const
    {
        uint64_2_t hi = {
            (a.hi & (~this->mask)) >> this->n,
            (a.lo & (~this->mask)) >> this->n
        };

        uint64_2_t lo = {
            a.hi & this->mask,
            a.lo & this->mask
        };

        uint64_2_t r = this->mul(hi, this->q_plus);
        r = { r.hi >> this->n, r.lo >> this->n };
        r = this->mul(r, this->mod_ast);
        r = { r.hi & this->mask, r.lo & this->mask };
        return this->add(r, lo);
    }

    uint64_2_t mont_mul(uint64_2_t a, uint64_2_t b) const
    {
        uint64_2_t t = this->mul(a, b);
        uint64_2_t u = this->mul(t, this->n_prime);

        u.hi &= this->mask;
        u.lo &= this->mask;
        uint64_2_t c = this->add(t, this->mul(u, {0, this->mod}));

        c.hi >>= this->n;
        c.lo >>= this->n;
        return c;
    }

    uint64_2_t add(uint64_2_t a, uint64_2_t b) const
    {
        uint64_t carry = a.lo & b.lo;
        return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
    }

    uint64_2_t negate(uint64_2_t a) const
    {
        return {
            a.lo ^ a.hi,
            a.lo
        };
    }

    uint64_2_t mul(uint64_2_t a, uint64_2_t b) const
    {
        return fast_mul(a, b);
    }

    uint64_2_t ref_mul(uint64_2_t a, uint64_2_t b) const
    {
        uint64_2_t c = { 0, 0 };

        for (int i = 0; i <= global::E.get_n(); i++)
        {
            uint64_2_t tmp = { 0, 0 };
            if ((a.hi >> i) & 1)
                tmp.hi = this->mask;
            if ((a.lo >> i) & 1)
                tmp.lo = this->mask;

            /* 2 bit carryless multiplier */
            uint64_2_t aib = {
                (b.lo & tmp.hi) ^ (b.hi & tmp.lo),
                b.lo & tmp.lo
            };

            aib.lo <<= i;
            aib.hi <<= i;
            c = global::E.add(c, aib);
        }

        return c;
    }

    uint64_2_t fast_mul(uint64_2_t a, uint64_2_t b) const
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
        for (int i = 0; i <= global::E.get_n(); i++)
        {
            if ((b.lo >> i)&1)
            {
                hi ^= (a.lo << i) & lo;
                lo ^= (a.lo << i);
            }
        }

        return { hi1 ^ hi2 ^ hi, lo };
    }

    int get_n() const { return this->n; }
    uint64_t get_mod() const { return this->mod; }
    uint64_t get_mask() const { return this->mask; }
};


class Extension_element
{
private:
    uint64_2_t repr;

public:
    Extension_element() { };

    Extension_element(const uint64_t lo, const uint64_t hi)
    {
        this->repr = { hi, lo };
    }

    Extension_element(const uint64_2_t repr)
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
        /* turn other to the additive inverse and then just add */
        return Extension_element(
            global::E.add(this->repr, global::E.negate(other.get_repr()))
        );
    }

    Extension_element &operator-=(const Extension_element &other)
    {
        this->repr = global::E.add(this->repr, global::E.negate(other.get_repr()));
        return *this;
    }

    Extension_element operator*(const Extension_element &other) const
    {
        uint64_2_t prod = global::E.mul(this->repr, other.get_repr());
        /* use rem in initializer? same for GF */
        return Extension_element(global::E.rem(prod));
    }

    Extension_element &operator*=(const Extension_element &other)
    {
        uint64_2_t prod = global::E.mul(this->repr, other.get_repr());
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
    uint64_2_t get_repr() const { return this->repr; }

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

#endif
