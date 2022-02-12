#include <stdint.h>
#include <immintrin.h>

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

bool GF2n::operator==(const GF2n &other) const
{
    return this->n == other.n && this->mod == other.mod;
}

GF_element GF2n::zero()
{
    return GF_element(0, *this);
}

GF_element GF2n::one()
{
    return GF_element(1, *this);
}

/* this can create zero, is it a problem? */
GF_element GF2n::random()
{
    return GF_element(global::randgen() & this->mask, *this);
}


/* GF element */

GF_element::GF_element(int64_t n, const GF2n &field) : field(field),
                                                       repr(n)
{
    return;
}

GF_element GF_element::operator+(const GF_element &other)
{
    return GF_element(this->repr ^ other.repr, this->field);
}

GF_element GF_element::operator*(const GF_element &other)
{
    const __m128i prod = _mm_clmulepi64_si128(
        _mm_set_epi64x(0, this->repr),
        _mm_set_epi64x(0, other.repr),
        0x0
    );

    uint64_t lo = _mm_extract_epi64(prod, 0x0);
    /* discard hi, only support up to 32 bit */

    return GF_element(
        util::modz2(lo, this->field.mod, this->field.n),
        this->field
    );
}

bool GF_element::operator==(const GF_element &other)
{
    return this->repr == other.repr && this->field == other.field;
}
