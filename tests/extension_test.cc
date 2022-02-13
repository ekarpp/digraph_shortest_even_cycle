#include <iostream>

#include "extension_test.hh"
#include "../src/extension.hh"
#include "../src/util.hh"
#include "../src/global.hh"

Extension global::E;
using namespace std;


Extension_test::Extension_test(int deg)
{
    this->n = deg;
    global::E.init(this->n, util::irred_poly(this->n));
    cout << "-----------------" << endl;
    cout << "TESTING EXTENSION" << endl;
    cout << "-----------------" << endl << endl;
}

void Extension_test::end_test(int err)
{
    if (err)
        cout << err << " errors";
    else
        cout << "CLEAR";
    cout << endl << endl;
}

void Extension_test::test_add_inverse()
{
    cout << "add inverse" << endl;
    int err = 0;
    for (int i = 0; i < 10000; i++)
    {
        Extension_element e = global::E.random();
        if (e - e != global::E.zero())
            err++;
    }
    end_test(err);
}
