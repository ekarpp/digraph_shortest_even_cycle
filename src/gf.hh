#ifndef GF_H
#define GF_H

#include <stdint.h>

class GF_element;

/* GF(2^n) */
class GF2n
{
private:
    /* largest possible element in the field */
    int64_t mask;
    int n;
    int64_t mod;

public:
    GF2n() {};
    void init(const int n, const int64_t mod);
    bool operator==(const GF2n &other) const;
    GF_element zero();
    GF_element one();
    GF_element random();
    int64_t rem(int64_t a) const;

    int get_n() const { return this->n; }
    int64_t get_mod() const { return this->mod; }
    int64_t get_mask() const { return this->mask; }
};

class GF_element
{
private:
    const GF2n &field;
    int64_t repr;

public:
    GF_element(const int64_t n, const GF2n &field);
    GF_element operator+(const GF_element &other);
    GF_element operator*(const GF_element &other);
    bool operator==(const GF_element &other);

    const GF2n &get_field() const { return this->field; }
    int64_t get_repr() const { return this->repr; }

    GF_element operator-(const GF_element &other)
    {
        return *this + other;
    }

    GF_element &operator=(const GF_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }
};



#endif
