#include "global.hh"
#include "extension.hh"


/* Extension */

void Extension::init(const int n, const int64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->mask = (1 << this->n) - 1;
    return;
}

Extension_element Extension::zero()
{
    return Extension_element(0b0, 0b0);
}

Extension_element Extension::one() const
{
    return Extension_element(0b1, 0b0);
}

Extension_element Extension::random()
{
    return Extension_element(
        global::randgen() & this->mask,
        global::randgen() & this->mask
    );
}

inline int64_2_t Extension::add(int64_2_t a, int64_2_t b) const
{
    int64_t carry = a.lo & b.lo;
    return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
}

/* mul if need 64b support, maybe store intermadiary
results in  array of len 64 and then add with offset */

/* Extension element */
Extension_element::Extension_element(const int64_t lo, const int64_t hi)
{
    this->repr = { hi, lo };
    return;
}

Extension_element::Extension_element(const int64_2_t repr)
{
    this->repr = repr;
    return;
}

Extension_element Extension_element::operator+(const Extension_element &other)
{
    return Extension_element(
        global::E.add(this->repr, other.get_repr())
    );
}

Extension_element Extension_element::operator-(const Extension_element &other)
{
    /* turn other to the additive inverse and then just add */
    int64_2_t b = {
        other.get_lo() ^ other.get_hi(),
        other.get_lo()
    };
    return Extension_element(
        global::E.add(this->repr, b)
    );
}

Extension_element Extension_element::operator*(const Extension_element &other)
{
    return global::E.one();
}

bool Extension_element::operator==(const Extension_element &other)
{
    return this->repr.lo == other.get_lo() && this->repr.hi == other.get_hi();
}
