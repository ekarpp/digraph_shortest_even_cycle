#include <iostream>
#include <vector>
#include <valarray>

#include "global.hh"
#include "graph.hh"
#include "gf.hh"

using namespace std;

Graph::Graph(vector<vector<int>> adjacency_list)
{
    n = adjacency_list.size();
    this->sample_adjacency(adjacency_list);

    cout << "created graph of " << n << " vertices:" << endl;
    for (uint i = 0; i < adjacency_list.size(); i++)
    {
        cout << i << ": ";
        for (uint j = 0; j < adjacency_list[i].size(); j++)
            cout << adjacency_list[i][j] << " ";
        cout << endl;
    }
}

/* samples the adjacency matrix with random edge weights from F
 * also creates a loop at each vertex
 */
void Graph::sample_adjacency(vector<vector<int>> adjacency_list)
{
    valarray<GF_element> m(global::F.zero(), this->n * this->n);

    for (int u = 0; u < this->n; u++)
    {
        /* loop at each vertex */
        m[u*this->n + u] = global::F.random();
        for (uint i = 0; i < adjacency_list[u].size(); i++)
        {
            int v = adjacency_list[u][i];
            m[u*this->n + v] = global::F.random();
        }
    }

    this->A = FMatrix(this->n, m);
    return;
}
