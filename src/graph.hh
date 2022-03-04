/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GRAPH_H
#define GRAPH_H
#include <vector>

#include "gf.hh"
#include "fmatrix.hh"

class Graph
{
private:
    int n;
    std::vector<std::vector<int>> adj;
    FMatrix A;

    void sample_adjacency();

public:
    Graph(std::vector<std::vector<int>> &adjacency_list);

    int get_n() const { return n; }
    FMatrix &get_A() { return A; }

    void dfs_cycle(int start, int *len) const;
};
#endif
