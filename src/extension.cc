#include "global.hh"
#include "extension.hh"


/* Extension */

void Extension::init(const int n, const int64_t mod)
{
    this->n = n;
    this->mod = mod;
    this->masklo = (1 << this->n) - 1;
    this->maskhi = this->masklo << 32;
    this->mask = this->masklo | this->maskhi;
    return;
}

Extension_element Extension::zero()
{
    return Extension_element(0);
}

Extension_element Extension::one() const
{
    return Extension_element(1);
}

Extension_element Extension::random()
{
    return Extension_element(global::randgen() & this->mask);
}

int64_t Extension::add(int64_t a, int64_t b) const
{
    int64_t result = a ^ b;
    int64_t carry = (a & b & this->masklo) << 32;
    return result ^ carry;
}

/* Extension element */
Extension_element::Extension_element(const int64_t n)
{
    this->repr = n;
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
    int64_t b = (other.get_lo() ^ other.get_hi()) << 32;
    b |= other.get_lo();
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
    return this->repr == other.repr;
}

int64_t Extension_element::get_lo() const
{
    return this->repr & global::E.get_masklo();
}

int64_t Extension_element::get_hi() const
{
    return (this->repr & global::E.get_maskhi()) >> 32;
}
