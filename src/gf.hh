#ifndef GF_H
#define GF_H

#include <stdint.h>

#include "algebra.hh"
#include "element.hh"

class GF2n;

class GF_element
{
private:
    const GF2n &field;

public:
    const int64_t repr;
    GF_element(const int64_t n, const GF2n &field);
    GF_element operator+(const GF_element &other);
    GF_element operator*(const GF_element &other);
    bool operator==(const GF_element &other);

    GF_element operator-(const GF_element &other)
    {
        return *this + other;
    }
};

/* GF(2^n) */
class GF2n
{
private:
    int n;
    int64_t mod;
public:
    GF2n() {};
    void init(const int n, const int64_t mod);
    bool operator==(const GF2n &other) const;
    GF_element zero();
    GF_element one();
    GF_element random();
};

#endif
