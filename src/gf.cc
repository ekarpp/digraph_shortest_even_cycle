/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <immintrin.h>
#include <iostream>

#include "extension.hh"
#include "gf.hh"
#include "global.hh"

/* GF */

void GF2n::init(const int n, const uint64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1ll << this->n) - 1;

    // 4.2 in https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011
    this->q_plus = this->quo(1ull << (2*this->n), mod);
    this->mod_ast = this->mask & mod;

    cout << "initialized GF(2^" << n << ") with modulus: ";
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

GF_element GF2n::zero() const
{
    return GF_element(0);
}

GF_element GF2n::one() const
{
    return GF_element(1);
}

/* this can create zero, is it a problem? */
GF_element GF2n::random() const
{
    return GF_element(global::randgen() & this->mask);
}

/* returns r s.t. for some q,
 * a = q*field.mod + r is the division relation (in Z(2^n))
 */
// a needs to be 128 bit for support up to GF(2^64)
// now just GF(2^32)
// 4.2 in https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011
uint64_t GF2n::rem(uint64_t a) const
{
    uint64_t lo = a & this->mask;
    uint64_t hi = (a & (~this->mask)) >> this->n;

    uint64_t rem = this->clmul(hi, this->q_plus);
    rem &= ~this->mask;
    rem >>= this->n;
    rem = this->clmul(rem, this->mod_ast);
    rem &= this->mask;
    rem ^= lo;
    return rem;
}

/* returns q s.t. for some r,
 * a = q*b + r is the division relation
 */
uint64_t GF2n::quo(uint64_t a, uint64_t b) const
{
    uint64_t q = 0b0;
    int degb;
    for (degb = 0; b >> degb; degb++);
    degb--;
    while (a >= b)
    {
        int shift;
        for (shift = 0; a >> shift; shift++);
        shift -= 1 + degb;
        /* shift = deg(a) - deg(b) */
        q ^= (1ll << shift);
        a ^= (b << shift);
    }
    return q;
}

/* carryless multiplication of a and b, polynomial multiplicatoin that is
 * done with Intel CLMUL
 */
uint64_t GF2n::clmul(uint64_t a, uint64_t b) const
{
    const __m128i prod = _mm_clmulepi64_si128(
        _mm_set_epi64x(0, a),
        _mm_set_epi64x(0, b),
        0x0
    );

    uint64_t lo = _mm_extract_epi64(prod, 0x0);
    /* discard hi, only support up to 32 bit */
    return lo;
}

/* GF element */

GF_element::GF_element(const uint64_t n)
{
    this->repr = n;
    return;
}

GF_element::GF_element(const GF_element &e)
{
    this->repr = e.get_repr();
}

GF_element GF_element::operator+(const GF_element &other) const
{
    return GF_element(this->repr ^ other.get_repr());
}

GF_element &GF_element::operator+=(const GF_element &other)
{
    this->repr ^= other.get_repr();
    return *this;
}

GF_element GF_element::operator*(const GF_element &other) const
{
    const uint64_t prod = global::F.clmul(
        this->repr,
        other.get_repr()
    );

    return GF_element(
        global::F.rem(prod)
    );
}

GF_element &GF_element::operator*=(const GF_element &other)
{
    const uint64_t prod = global::F.clmul(
        this->repr,
        other.get_repr()
    );

    this->repr = global::F.rem(prod);

    return *this;
}


/* invert elements with Fermat's little theorem
 * and binary exponentiation.
 * From FLT we have: a^{p-1} = 1 mod p
 * <=> a^{p-2} = a^{-1} mod p.
 * with binary exponentiation:
 * a^{p-2} = a^2 a^{2^2} ... a^{2^{p-1}} */
GF_element GF_element::inv() const
{
    GF_element inv = GF_element(*this);

    for (int i = 0; i < global::F.get_n() - 2; i++)
    {
        inv *= inv;
        inv *= *this;
    }
    inv *= inv;
    return inv;
}

GF_element GF_element::operator/(const GF_element &other) const
{
    return *this * other.inv();
}

GF_element &GF_element::operator/=(const GF_element &other)
{
    GF_element inv = other.inv();
    const uint64_t prod = global::F.clmul(
        this->repr,
        inv.get_repr()
    );

    this->repr = global::F.rem(prod);
    return *this;
}

bool GF_element::operator==(const GF_element &other) const
{
    return this->repr == other.get_repr();
}

Extension_element GF_element::lift() const
{
    return Extension_element(this->repr, 0b0);
}
