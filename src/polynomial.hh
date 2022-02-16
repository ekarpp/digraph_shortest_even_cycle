/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>

#include "gf.hh"

class Polynomial
{
private:
    std::vector<GF_element> coeffs;
    int deg;

public:
    Polynomial(int n);

    Polynomial &operator*=(GF_element val);

    Polynomial &operator+=(const Polynomial &other);

    GF_element operator[](int i) const
    {
        return this->coeffs[i];
    }

    /* set coefficient with deg i to val */
    void operator()(int i, GF_element val)
    {
        this->coeffs[i] = val;
    }
};

#endif
