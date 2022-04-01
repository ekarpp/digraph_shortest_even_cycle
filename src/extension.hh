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

    uint64_2_t fast_mul(uint64_2_t a, uint64_2_t b) const;

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
