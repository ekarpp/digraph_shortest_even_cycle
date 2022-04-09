/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "gf.hh"
#include "global.hh"
#include "extension.hh"
#include "util.hh"

using namespace std;

/* Extension */
void Extension::init(const int n, const uint64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1ll << this->n) - 1;

    this->mod_ast = { 0, this->mod & this->mask };
    this->q_plus = this->quo({0, 1ull << (2*this->n)} , { 0, this->mod });

    cout << "initialized E(4^" << n << ") with modulus: ";
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

Extension_element Extension::zero() const
{
    return Extension_element(0b0, 0b0);
}

Extension_element Extension::one() const
{
    return Extension_element(0b1, 0b0);
}

Extension_element Extension::random() const
{
    return Extension_element(
        global::randgen() & this->mask,
        global::randgen() & this->mask
    );
}

/* Extension element */
GF_element Extension_element::project() const
{
    return GF_element(this->repr.lo);
}
