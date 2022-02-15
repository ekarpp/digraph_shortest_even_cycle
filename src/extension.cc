#include <iostream>

#include "gf.hh"
#include "global.hh"
#include "extension.hh"


/* Extension */

void Extension::init(const int n, const uint64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1 << this->n) - 1;

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

uint64_2_t Extension::rem(uint64_2_t a) const
{
    while (a.lo > this->mask || a.hi > this->mask)
    {
        int shift;
        for (shift = 0; a.lo >> shift || a.hi >> shift; shift++);
        shift -= 1 + this->n;
        /* mod has coefficients 01, thus its negation
         * is just it applied to hi and lo (see negate function) */
        a = this->add(a, { this->mod << shift, this->mod << shift });
    }
    return a;
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

/* multiplication by constant, 0 <= a < 4 */
uint64_2_t Extension::mul_const(int a, uint64_2_t b) const
{
    uint64_2_t c = { 0, 0 };
    /* form repeating polynomial and then
     * multiply each coefficient */
    if (a & 0b1)
        c.lo = this->mask;
    if (a & 0b10)
        c.hi = this->mask;

    /* lo bit if multiplication is simply and */
    uint64_t lo = c.lo & b.lo;
    /* boolean formula for the hi bit for multiplication
     * of two mod 4 integers. solved by writing the truth
     * table and then forming DNF and simplifying it
     * formula for a*b: (hi_a & lo_b & ~hi_b) |
     * (hi_a & lo_b & ~lo_a) | (lo_a & hi_b & ~lo_b) |
     * (lo_a & hi_b & ~hi_a)
     */
    uint64_t hi = c.hi & b.lo & (b.hi ^ this->mask);
    hi |= c.hi & b.lo & (c.lo ^ this->mask);
    hi |= c.lo & b.hi & (b.lo ^ this->mask);
    hi |= c.lo & b.hi & (c.hi ^ this->mask);

    return { hi, lo };
}

uint64_2_t Extension::mul(uint64_2_t a, uint64_2_t b) const
{
    /* this is horrible all around
     * how to make better?
     * (optimiza boolean formula with XOR?,
     * lookuptable for repeating constants?) */
    uint64_2_t c = { 0, 0 };

    for (int i = 0; i <= global::E.get_n(); i++)
    {
        int hi = (a.hi >> i) & 1;
        int lo = (a.lo >> i) & 1;
        uint64_2_t aib = global::E.mul_const((hi << 1) | lo, b);
        aib.lo <<= i;
        aib.hi <<= i;
        c = global::E.add(c, aib);
    }

    return c;
}

/* Extension element */
Extension_element::Extension_element(const uint64_t lo, const uint64_t hi)
{
    this->repr = { hi, lo };
    return;
}

Extension_element::Extension_element(const Extension_element &e)
{
    this->repr = { e.get_hi(), e.get_lo() };
}

Extension_element::Extension_element(const uint64_2_t repr)
{
    this->repr = repr;
    return;
}

Extension_element Extension_element::operator+(const Extension_element &other) const
{
    return Extension_element(
        global::E.add(this->repr, other.get_repr())
    );
}

Extension_element &Extension_element::operator+=(const Extension_element &other)
{
    this->repr = global::E.add(this->repr, other.get_repr());
    return *this;
}

Extension_element Extension_element::operator-(const Extension_element &other) const
{
    /* turn other to the additive inverse and then just add */
    return Extension_element(
        global::E.add(this->repr, global::E.negate(other.get_repr()))
    );
}

Extension_element &Extension_element::operator-=(const Extension_element &other)
{
    this->repr = global::E.add(this->repr, global::E.negate(other.get_repr()));
    return *this;
}

GF_element Extension_element::project() const
{
    return GF_element(this->repr.lo);
}

Extension_element Extension_element::operator*(const Extension_element &other) const
{
    uint64_2_t prod = global::E.mul(this->repr, other.get_repr());
    /* use rem in initializer? same for GF */
    return Extension_element(global::E.rem(prod));
}

Extension_element &Extension_element::operator*=(const Extension_element &other)
{
    uint64_2_t prod = global::E.mul(this->repr, other.get_repr());
    this->repr = global::E.rem(prod);
    return *this;
}

bool Extension_element::operator==(const Extension_element &other) const
{
    return this->repr.lo == other.get_lo() && this->repr.hi == other.get_hi();
}
