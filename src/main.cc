/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <sstream>

#include "gf.hh"
#include "global.hh"
#include "graph.hh"
#include "util.hh"
#include "fmatrix.hh"
#include "solver.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;

void parse_file(string fname, Graph &G)
{
    string line;
    ifstream file(fname);

    if (!file.is_open())
    {
        cout << "unable to open file: " << fname << endl;
        return;
    }

    vector<vector<int>> graph;
    int u = 0;
    while (getline(file, line))
    {
        /* comment */
        if (line[0] == '#')
            continue;
        istringstream iss(line);
        int v;
        vector<int> vec;
        while (iss >> v)
            vec.push_back(v);
        graph.push_back(vec);
        u++;
    }
    file.close();

    G.init(graph);
}

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "provide path to graph file with -f" << endl;
        cout << "line i (starting at zero) in graph file" << endl;
        cout << "tells to which nodes there is an edge to" << endl;
        cout << "-b to use the brute force solver instead (SLOW)" << endl;
        cout << "-q for no progress output from solver" << endl;
        return 0;
    }

    int opt;
    Graph G;

    bool brute = false;
    bool output = true;

    while ((opt = getopt(argc, argv, "qbf:")) != -1)
    {
        switch (opt)
        {
        case 'f':
            parse_file(optarg, G);
            /* error during parsing */
            if (G.get_n() == -1)
                return -1;
            break;
        case 'b':
            brute = true;
            break;
        case 'q':
            output = false;
            break;
        case '?':
            cout << "call with no arguments for help" << endl;
            return -1;
        }
    }

    uint64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    Solver s;

    if (brute)
        cout << s.shortest_even_cycle_brute(G) << endl;
    else
        cout << s.shortest_even_cycle(G, output) << endl;

    return 0;
}
