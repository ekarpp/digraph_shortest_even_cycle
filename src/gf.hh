/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GF_H
#define GF_H

#include <stdint.h>
#include <bitset>
#include <iostream>

#include "extension.hh"

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
    uint64_t quo(uint64_t a, uint64_t b) const;

    uint64_t q_plus;
    uint64_t mod_ast;

public:
    GF2n() {};
    void init(const int n, const uint64_t mod);
    GF_element zero() const;
    GF_element one() const;
    GF_element random() const;
    uint64_t rem(uint64_t a) const;
    uint64_t clmul(uint64_t a, uint64_t) const;

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
    GF_element(const uint64_t n);
    GF_element(const GF_element &e);
    GF_element operator+(const GF_element &other) const;
    GF_element &operator+=(const GF_element &other);
    GF_element operator*(const GF_element &other) const;
    GF_element &operator*=(const GF_element &other);
    GF_element operator/(const GF_element &other) const;
    GF_element &operator/=(const GF_element &other);
    bool operator==(const GF_element &other) const;

    GF_element inv() const;
    Extension_element lift() const;

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
