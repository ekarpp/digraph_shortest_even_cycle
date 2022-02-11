#include "graph.hh"
#include <vector>
#include <iostream>

using namespace std;

Graph::Graph(vector<vector<int>> adjacency_list)
{
    n = adjacency_list.size();
    sample_adjacency(adjacency_list);

    cout << "created graph of " << n << " vertices:" << endl;
    for (uint i = 0; i < adjacency_list.size(); i++)
    {
        cout << i << ": ";
        for (uint j = 0; j < adjacency_list[i].size(); j++)
            cout << adjacency_list[i][j] << " ";
        cout << endl;
    }
}

/* samples the adjacency matrix with random edge weights */
/* also creates a loop at each vertex */
void Graph::sample_adjacency(vector<vector<int>> adjacency_list)
{
    return;
}
