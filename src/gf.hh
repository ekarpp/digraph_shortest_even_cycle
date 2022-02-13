#ifndef GF_H
#define GF_H

#include <stdint.h>

#include "extension.hh"

/* forward declare */
class GF_element;
class Extension_element;

/* GF(2^n) */
class GF2n
{
private:
    /* largest possible element in the field */
    int64_t mask;
    int n;
    int64_t mod;
    int64_t quo(int64_t a, int64_t b) const;

public:
    GF2n() {};
    void init(const int n, const int64_t mod);
    GF_element zero() const;
    GF_element one() const;
    GF_element random();
    int64_t rem(int64_t a) const;
    int64_t ext_euclid(int64_t a) const;
    int64_t clmul(int64_t a, int64_t) const;

    int get_n() const { return this->n; }
    int64_t get_mod() const { return this->mod; }
    int64_t get_mask() const { return this->mask; }
};

class GF_element
{
private:
    int64_t repr;

public:
    GF_element(const int64_t n);
    GF_element operator+(const GF_element &other);
    GF_element operator*(const GF_element &other);
    GF_element operator/(const GF_element &other);
    bool operator==(const GF_element &other);

    Extension_element lift();

    int64_t get_repr() const { return this->repr; }

    GF_element operator-(const GF_element &other)
    {
        return *this + other;
    }

    GF_element &operator=(const GF_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    bool operator!=(const GF_element &other)
    {
        return !(*this == other);
    }

};

#endif
