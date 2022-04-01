/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "gf.hh"
#include "global.hh"
#include "extension.hh"
#include "util.hh"

using namespace std;

/* Extension */
void Extension::init(const int n, const uint64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1ll << this->n) - 1;

    this->mod_ast = { 0, this->mod & this->mask };
    this->q_plus = this->quo({0, 1ull << (2*this->n)} , { 0, this->mod });

    cout << "initialized E(4^" << n << ") with modulus: ";
    for (int i = n; i >= 0; i--)
    {
        if ((mod >> i) & 1)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;

    return;
}

Extension_element Extension::zero() const
{
    return Extension_element(0b0, 0b0);
}

Extension_element Extension::one() const
{
    return Extension_element(0b1, 0b0);
}

Extension_element Extension::random() const
{
    return Extension_element(
        global::randgen() & this->mask,
        global::randgen() & this->mask
    );
}

/* euclidean division, only used once during initialization.
 * b has to be monic for this to work */
uint64_2_t Extension::quo(uint64_2_t a, uint64_2_t b) const
{
    uint64_2_t q = { 0, 0 };
    int degb = this->n;
    int dega = 63 - min(__builtin_clzl(a.lo), __builtin_clzl(a.hi));
    while (dega >= degb)
    {
        uint64_2_t s = {
            (a.hi & (1ll << dega)) >> degb,
            (a.lo & (1ll << dega)) >> degb
        };

        q = this->add(q, s);
        a = this->add(a, this->negate(this->mul(b, s)));

        dega = 63 - min(__builtin_clzl(a.lo), __builtin_clzl(a.hi));
    }

    return q;
}

/* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
uint64_2_t Extension::rem(uint64_2_t a) const
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

uint64_2_t Extension::add(uint64_2_t a, uint64_2_t b) const
{
    uint64_t carry = a.lo & b.lo;
    return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
}

uint64_2_t Extension::negate(uint64_2_t a) const
{
    return {
        a.lo ^ a.hi,
        a.lo
    };
}

uint64_2_t Extension::mul(uint64_2_t a, uint64_2_t b) const
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

/* Extension element */
GF_element Extension_element::project() const
{
    return GF_element(this->repr.lo);
}
