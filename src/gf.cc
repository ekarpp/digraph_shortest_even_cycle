/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <iostream>

#include "extension.hh"
#include "gf.hh"
#include "global.hh"
#include "util.hh"

using namespace std;

/* GF */
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
