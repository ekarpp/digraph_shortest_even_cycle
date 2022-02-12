#include <iostream>

#include "gf_test.hh"
#include "../src/gf.hh"
#include "../src/util.hh"

using namespace std;

GF_test::GF_test(int deg)
{
    GF2n field();
    this->n = deg;
    this->field.init(this->n, util::irred_poly(this->n));
}

void GF_test::end_test(int err)
{
    if (err)
        cout << err << " errors";
    else
        cout << "CLEAR";
    cout << endl << endl;
}

void GF_test::test_add_inverse()
{
    cout << "add inverse" << endl;
    int err = 0;
    int64_t i = 0;
    while (i <= this->field.get_mask())
    {
        GF_element e(i, this->field);
        if (e + e != this->field.zero()
            || e - e != this->field.zero())
            err++;
        i++;
    }
    end_test(err);
}

void GF_test::test_associativity()
{
    cout << "test associativity" << endl;
    int err = 0;
    for (int i = 0; i < 10000; i++)
    {
        GF_element a = this->field.random();
        GF_element b = this->field.random();
        GF_element c = this->field.random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    this->end_test(err);
}

void GF_test::test_mul_id()
{
    cout << "mul with id" << endl;
    int err = 0;
    int64_t i = 0;
    while (i <= this->field.get_mask())
    {
        GF_element e(i, this->field);
        if (e * this->field.one() != e)
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_mul_inverse()
{
    cout << "mul with inverse" << endl;
    int err = 0;
    int64_t i = 2;
    while (i <= this->field.get_mask())
    {
        GF_element e(i, this->field);
        if (e / e != this->field.one())
            err++;
        i++;
    }
    this->end_test(err);
}
