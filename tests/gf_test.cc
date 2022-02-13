#include <iostream>

#include "gf_test.hh"
#include "../src/gf.hh"
#include "../src/util.hh"
#include "../src/global.hh"
#include "../src/extension.hh"

using namespace std;
GF2n global::F;

GF_test::GF_test(int deg)
{
    this->n = deg;
    global::F.init(this->n, util::irred_poly(this->n));
    cout << "----------" << endl;
    cout << "TESTING GF" << endl;
    cout << "----------" << endl << endl;
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
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e + e != global::F.zero()
            || e - e != global::F.zero())
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
        GF_element a = global::F.random();
        GF_element b = global::F.random();
        GF_element c = global::F.random();
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
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e * global::F.one() != e)
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_mul_inverse()
{
    cout << "mul with inverse" << endl;
    int err = 0;
    int64_t i = 1;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e / e != global::F.one())
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_lift_project()
{
    cout << "lift project" << endl;
    int err = 0;
    int64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        Extension_element b(i, global::randgen() & global::E.get_mask());
        Extension_element c(i, 0b0);
        if (e.lift().project() != e
            || b.project() != e
            || e.lift() != c)
            err++;
        i++;
    }
    this->end_test(err);
}
