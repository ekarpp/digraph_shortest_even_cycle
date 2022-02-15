/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_H
#define EXTENSION_H

#include <stdint.h>

#include "gf.hh"

/* forward declare */
class GF_element;
class Extension_element;

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
/* vectorize with AVX? */
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
    /* lookuptable that returns for each non-negative integer < 4
     * the polynomial where each coefficient is that integer */
    uint64_2_t lookup[4];

public:
    Extension() {};
    void init(const int n, const uint64_t mod);
    Extension_element zero() const;
    Extension_element one() const;
    Extension_element random() const;

    uint64_2_t rem(uint64_2_t a) const;
    uint64_2_t add(uint64_2_t a, uint64_2_t b) const;
    uint64_2_t mul_const(int a, uint64_2_t b) const;
    uint64_2_t mul(uint64_2_t a, uint64_2_t b) const;
    uint64_2_t negate(uint64_2_t a) const;

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
    Extension_element(const uint64_t lo, const uint64_t hi);
    Extension_element(const uint64_2_t repr);
    Extension_element(const Extension_element& e);
    Extension_element operator+(const Extension_element &other) const;
    Extension_element &operator+=(const Extension_element &other);
    Extension_element operator-(const Extension_element &other) const;
    Extension_element &operator-=(const Extension_element &other);
    Extension_element operator*(const Extension_element &other) const;
    Extension_element &operator*=(const Extension_element &other);
    bool operator==(const Extension_element &other) const;

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
};

#endif
