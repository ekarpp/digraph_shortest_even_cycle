/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GRAPH_H
#define GRAPH_H
#include <vector>

#include "gf.hh"
#include "fmatrix.hh"

using namespace std;

class Graph
{
private:
    int n;
    void sample_adjacency(vector<vector<int>> adjacency_list);
    FMatrix A;

public:
    Graph(vector<vector<int>> adjacency_list);

    int get_n() const { return n; }
    FMatrix &get_A() { return A; }
};
#endif
