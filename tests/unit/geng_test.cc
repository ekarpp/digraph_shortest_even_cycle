/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>

#include "geng_test.hh"
#include "../../src/graph.hh"

using namespace std;

Geng_test::Geng_test()
{
    cout << "------------" << endl;
    cout << "TESTING GENG" << endl;
    cout << "------------" << endl;

    this->run();
}

void Geng_test::test_geng()
{
    string line = "";

    while (cin >> line)
    {
        cout << line << endl;
    }
}
