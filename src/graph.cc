/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <valarray>
#include <stack>

#include "global.hh"
#include "graph.hh"
#include "gf.hh"

using namespace std;

Graph::Graph(vector<vector<int>> &adjacency_list)
{
    this->n = adjacency_list.size();
    this->adj = adjacency_list;
    this->sample_adjacency();

    cout << "created graph of " << this->n << " vertices:" << endl;
    /* write this to file ?? */
    /*
    for (uint i = 0; i < adjacency_list.size(); i++)
    {
        cout << i << ": ";
        for (uint j = 0; j < adjacency_list[i].size(); j++)
            cout << adjacency_list[i][j] << " ";
        cout << endl;
    }
    */
}

/* samples the adjacency matrix with random edge weights from F
 * also creates a loop at each vertex
 */
void Graph::sample_adjacency()
{
    valarray<GF_element> m(global::F.zero(), this->n * this->n);

    for (int u = 0; u < this->n; u++)
    {
        /* loop at each vertex */
        m[u*this->n + u] = global::F.random();
        for (uint i = 0; i < this->adj[u].size(); i++)
        {
            int v = this->adj[u][i];
            m[u*this->n + v] = global::F.random();
        }
    }

    this->A = FMatrix(this->n, m);
    return;
}

/* goes through all cycles that contain vertex start
 * and updates len accordingly.
 * len contains the length of the shortest found so far */
void Graph::dfs_cycle(int start, int *len) const
{
    int depth = 0;
    valarray<bool> visited(false, this->n);

    stack<int> S;
    S.push(start);

    while (!S.empty())
    {
        int v = S.top();
        S.pop();

        if (visited[v])
        {
            if (v == start && depth % 2 == 0 && depth < *len)
                *len = depth;
            depth--;
        }
        else
        {
            visited[v] = true;
            depth++;
            vector<int> nbors = this->adj[v];
            for (uint i = 0; i < nbors.size(); i++)
                S.push(nbors[i]);
        }
    }
}
