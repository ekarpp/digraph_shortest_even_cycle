#ifndef EXTENSION_H
#define EXTENSION_H

#include <stdint.h>

class Extension_element;

/* Extension of GF(2^n) to the ring E(4^n).
 * If GF(2^n) = Z2 / <g_2> for irreducible polynomial
 * g_2 of degree n, then if g_4 is g_2 but coefficients
 * projected to Z4 we have E(4^n) = Z4 / <g_4>. */
class Extension
{
private:
    int n;
    int64_t mod;

public:
    Extension() {};
    void init(const int n, const int64_t mod);
    Extension_element zero();
    Extension_element one();
    Extension_element random();

    int get_n() { return this->n; }
    int64_t get_mod() { return this->mod; }
};

class Extension_element
{
private:
    const Extension &ring;
    int64_t repr;

public:
    Extension_element(const int64_t n, const Extension &ring);
    Extension_element operator+(const Extension_element &other);
    Extension_element operator*(const Extension_element &other);
    Extension_element operator-(const Extension_element &other);
    bool operator==(const Extension_element &other);

    const Extension &get_ring() { return this->ring; }
    int64_t get_repr() { return this->repr; }

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    bool operator!=(const Extension_element &other)
    {
        return !(*this == other);
    }
}

#endif
