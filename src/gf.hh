/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GF_H
#define GF_H

#include <stdint.h>
#include <bitset>
#include <iostream>
#include <immintrin.h>

#include "extension.hh"
#include "global.hh"

/* forward declare */
class GF_element;
class Extension_element;

/* GF(2^n) */
class GF2n
{
private:
    /* largest possible element in the field */
    uint64_t mask;
    int n;
    uint64_t mod;

    /* returns q s.t. for some r,
     * a = q*b + r is the division relation
     */
    uint64_t quo(uint64_t a, uint64_t b) const
    {
        uint64_t q = 0b0;
        int degb = 63 - __builtin_clzl(b);
        while (a >= b)
        {
            int shift = 63 - __builtin_clzl(a);
            shift -= degb;
            /* shift = deg(a) - deg(b) */
            q ^= (1ll << shift);
            a ^= (b << shift);
        }
        return q;
    }

public:
    GF2n()
    {
        this->n = GF2_bits;

#if GF2_bits == 16
        /* x^16 + x^5 + x^3 + x^2 +  1 */
        this->mod = 0x1002D;
        this->mask = 0xFFFF;
#elif GF2_bits == 32
        /* x^32 + x^7 + x^3 + x^2 + 1 */
        this->mod = 0x10000008D;
        this->mask = 0xFFFFFFFF;
#else
        GF2_bits_eq_16_or_32
#endif
        if (global::output)
        {
            std::cout << "initialized GF(2^" << this->n << ") with modulus: ";
            for (int i = n; i >= 0; i--)
            {
                if ((this->mod >> i) & 1)
                    std::cout << "1";
                else
                    std::cout << "0";
            }
            std::cout << std::endl;
        }
    }

    GF_element zero() const;
    GF_element one() const;
    GF_element random() const;

    /* returns r s.t. for some q,
     * a = q*field.mod + r is the division relation (in Z(2^n))
     */
    uint64_t rem(uint64_t a) const
    {
        uint64_t lo = a & this->mask;
        uint64_t hi = a >> this->n;

#if GF2_bits == 16
        uint64_t rem = hi ^ (hi >> 14) ^ (hi >> 13) ^ (hi >> 11);
        rem ^= (rem << 2) ^ (rem << 3) ^ (rem << 5);
        rem ^= lo;
#elif GF2_bits == 32
        uint64_t rem = hi ^ (hi >> 30) ^ (hi >> 29) ^ (hi >> 25);
        rem ^= (rem << 2) ^ (rem << 3) ^ (rem << 7);
        rem ^= lo;
#endif
        return rem & this->mask;
    }

    /* returns s s.t. for some t: s*a + t*field.mod = gcd(field.mod, a)
     * <=> s*a + t*field.mod = 1 taking mod field.mod we get
     * s*a = 1 mod field.mod and thus a^-1 = s mod field.mod*/
    uint64_t ext_euclid(uint64_t a) const
    {
        // assert(a != 0)
        uint64_t s = 0x1;
        uint64_t s_next = 0x0;
        uint64_t r = a;
        uint64_t r_next = this->mod;
        uint64_t tmp;

        while (r_next != 0x0)
        {
            uint64_t q = this->quo(r, r_next);
            tmp = r ^ this->clmul(q, r_next);
            r = r_next;
            r_next = tmp;

            tmp = s ^ this->clmul(q, s_next);
            s = s_next;
            s_next = tmp;
        }

        return s;
    }

    /* carryless multiplication of a and b, polynomial multiplicatoin that is
     * done with Intel CLMUL
     */
    uint64_t clmul(uint64_t a, uint64_t b) const
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

    int get_n() const { return this->n; }
    uint64_t get_mod() const { return this->mod; }
    uint64_t get_mask() const { return this->mask; }
};

class GF_element
{
private:
    uint64_t repr;

public:
    GF_element() { }

    GF_element(const uint64_t n)
    {
        this->repr = n;
    }

    GF_element(const GF_element &e)
    {
        this->repr = e.get_repr();
    }

    GF_element operator+(const GF_element &other) const
    {
        return GF_element(this->repr ^ other.get_repr());
    }

    GF_element &operator+=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    GF_element operator*(const GF_element &other) const
    {
        const uint64_t prod = global::F.clmul(
            this->repr,
            other.get_repr()
        );

        return GF_element(
            global::F.rem(prod)
        );
    }

    GF_element &operator*=(const GF_element &other)
    {
        const uint64_t prod = global::F.clmul(
            this->repr,
            other.get_repr()
        );

        this->repr = global::F.rem(prod);

        return *this;
    }

    GF_element inv() const
    {
        return GF_element(global::F.ext_euclid(this->repr));
    }

    void inv_in_place()
    {
        this->repr = global::F.ext_euclid(this->repr);
    }

    GF_element operator/(const GF_element &other) const
    {
        return *this * other.inv();
    }

    GF_element &operator/=(const GF_element &other)
    {
        const uint64_t inv = global::F.ext_euclid(other.get_repr());
        const uint64_t prod = global::F.clmul(
            this->repr,
            inv
            );

        this->repr = global::F.rem(prod);
        return *this;
    }

    bool operator==(const GF_element &other) const
    {
        return this->repr == other.get_repr();
    }

    uint64_t get_repr() const { return this->repr; }

    GF_element operator-(const GF_element &other) const
    {
        return *this + other;
    }

    GF_element &operator-=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    GF_element &operator=(const GF_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    Extension_element lift() const;

    bool operator!=(const GF_element &other) const
    {
        return !(*this == other);
    }

    bool operator>(const GF_element &other) const
    {
        return this->repr > other.get_repr();
    }

    void print() const
    {
        std::cout << std::bitset<8>(this->repr) << std::endl;
    }
};

#endif
