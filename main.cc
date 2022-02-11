#include <iostream>
#include <vector>
#include "src/graph.hh"

using namespace std;

int main(int argc, char **argv)
{
    vector<vector<int>> vec = {
        {1},
        {0}
    };
    Graph G(vec);

    return 0;
}
