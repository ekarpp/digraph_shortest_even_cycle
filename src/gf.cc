/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <immintrin.h>
#include <iostream>

#include "extension.hh"
#include "gf.hh"
#include "global.hh"
#include "util.hh"

using namespace std;

/* GF */
void GF2n::init(const int n, const uint64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1ll << this->n) - 1;

    // 4.2 in https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011
    this->q_plus = this->quo(1ull << (2*this->n), mod);
    this->mod_ast = this->mask & mod;

    cout << "initialized GF(2^" << n << ") with modulus: ";
    for (int i = n; i >= 0; i--)
    {
        if ((mod >> i) & 1)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;

    return;
}

GF_element GF2n::zero() const
{
    return GF_element(0);
}

GF_element GF2n::one() const
{
    return GF_element(1);
}

/* this can create zero, is it a problem? */
GF_element GF2n::random() const
{
    return GF_element(global::randgen() & this->mask);
}

/* TODO: solve forward declaration issue and move to .hh */
Extension_element GF_element::lift() const
{
    return Extension_element(this->repr, 0b0);
}
