/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>

#include "geng_test.hh"
#include "../../src/graph.hh"
#include "../../src/solver.hh"

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
    cin >> line;
    int n = stoi(line);
    Solver s;

    int total = 0;
    int fail = 0;
    cout << "failed:" << endl;
    do
    {
        vector<vector<int>> graph;
        for (int i = 0; i < n; i++)
        {
            vector<int> vec;
            cin >> line;
            for (int j = 0; j < n; j++)
            {
                if (line[j] == '1')
                    vec.push_back(j);
            }
            graph.push_back(vec);

        }
        Graph G(graph, false);
        total++;
        int algo = s.shortest_even_cycle(G, false);
        int ref = s.shortest_even_cycle_brute(G);
        if (algo != ref)
            fail++;
        if (total % 1000 == 1)
        {
            cout << "\r" << fail << "/" << total;
            cout << " (" << ((float) fail) / total * 100 << "%)";
            cout << flush;
        }
    } while (cin >> line);
    cout << "\r" << fail << "/" << total;
    cout << " (" << ((float) fail) / total * 100 << "%)";
    cout << endl;
}
