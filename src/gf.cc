#include <stdint.h>
#include <immintrin.h>
#include <iostream>

#include "extension.hh"
#include "gf.hh"
#include "global.hh"

/* GF */

void GF2n::init(const int n, const int64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1 << this->n) - 1;
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
// r needs to be 128 bit for support up to GF(2^64)
// now just GF(2^32)
int64_t GF2n::rem(int64_t a) const
{
    while (a > this->mask)
    {
        int shift;
        for (shift = 0; a >> shift; shift++);
        shift -= 1 + this->n;
        /* shift = deg(a) - deg(b) */
        a ^= (this->mod << shift);
    }
    return a;
}

/* returns q s.t. for some r,
 * a = q*b + r is the division relation
 */
int64_t GF2n::quo(int64_t a, int64_t b) const
{
    int64_t q = 0b0;
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
int64_t GF2n::clmul(int64_t a, int64_t b) const
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

/* returns s s.t. for some t: s*a + t*field.mod = gcd(field.mod, a)
 * <=> s*a + t*field.mod = 1 taking mod field.mod we get
 * s*a = 1 mod field.mod and thus a^-1 = s mod field.mod*/
int64_t GF2n::ext_euclid(int64_t a) const
{
    // assert(a != 0)
    int64_t s = 0x1;
    int64_t s_next = 0x0;
    int64_t r = a;
    int64_t r_next = this->mod;
    int64_t tmp;

    while (r_next != 0x0)
    {
        int64_t q = this->quo(r, r_next);
        tmp = r ^ this->clmul(q, r_next);
        r = r_next;
        r_next = tmp;

        tmp = s ^ this->clmul(q, s_next);
        s = s_next;
        s_next = tmp;
    }

    return s;
}


/* GF element */

GF_element::GF_element(const int64_t n)
{
    this->repr = n;
    return;
}

GF_element GF_element::operator+(const GF_element &other) const
{
    return GF_element(this->repr ^ other.get_repr());
}

GF_element GF_element::operator*(const GF_element &other) const
{
    const int64_t prod = global::F.clmul(
        this->repr,
        other.get_repr()
    );

    return GF_element(
        global::F.rem(prod)
    );
}

GF_element GF_element::inv() const
{
    return GF_element(global::F.ext_euclid(this->repr));
}

GF_element GF_element::operator/(const GF_element &other) const
{
    return *this * other.inv();
}

bool GF_element::operator==(const GF_element &other) const
{
    return this->repr == other.get_repr();
}

Extension_element GF_element::lift() const
{
    return Extension_element(this->repr, 0b0);
}
