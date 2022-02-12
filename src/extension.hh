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
    int64_t masklo;
    int64_t maskhi;
    int64_t mask;

public:
    Extension() {};
    void init(const int n, const int64_t mod);
    Extension_element zero();
    Extension_element one() const;
    Extension_element random();

    int64_t add(int64_t a, int64_t b) const;

    int get_n() { return this->n; }
    int64_t get_mod() { return this->mod; }
    int64_t get_mask() { return this->mask; }
    int64_t get_masklo() const { return this->masklo; }
    int64_t get_maskhi() const { return this->maskhi; }
};

/* represented as aaaabbbb
 * where aaaa are high bits of each coefficient
 * and bbbb are the low bits
 */
class Extension_element
{
private:
    int64_t repr;

public:
    Extension_element(const int64_t n);
    Extension_element operator+(const Extension_element &other);
    Extension_element operator-(const Extension_element &other);
    Extension_element operator*(const Extension_element &other);
    bool operator==(const Extension_element &other);

    int64_t get_repr() const { return this->repr; }
    int64_t get_lo() const;
    int64_t get_hi() const;

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    bool operator!=(const Extension_element &other)
    {
        return !(*this == other);
    }
};

#endif
