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
        cout << s.shortest_even_cycle(G, false) << endl;
    } while (cin >> line);
}
