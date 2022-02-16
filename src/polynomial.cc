/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */

#include "polynomial.hh"
#include "global.hh"

/* be lazy and just store coefficients in vector of length n.
 * dont care if some of the coefficients are zero */
Polynomial::Polynomial(int n): coeffs(n+1, global::F.zero())
{
    this->deg = n;
}

Polynomial &Polynomial::operator*=(GF_element val)
{
    for (int i = 0; i < this->deg; i++)
        this->coeffs[i] *= val;

    return *this;
}

Polynomial &Polynomial::operator+=(const Polynomial &other)
{
    // assert(this->deg == other.deg)
    for (int i = 0; i < this->deg; i++)
        this->coeffs[i] += other[i];

    return *this;
}
