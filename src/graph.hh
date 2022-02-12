#ifndef GRAPH_H
#define GRAPH_H
#include <vector>

using namespace std;

class Graph
{
private:
    int n;
    void sample_adjacency(vector<vector<int>> adjacency_list);

public:
    Graph(vector<vector<int>> adjacency_list);

    int get_n() { return n; }
};
#endif
