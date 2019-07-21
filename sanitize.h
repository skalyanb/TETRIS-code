//
// Created by Suman Kalyan Bera on 2019-07-19.
//

#ifndef SUBGRAPHCOUNT_SANITIZE_H
#define SUBGRAPHCOUNT_SANITIZE_H

#include <cstdlib>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>

namespace sanitize
{
    using VertexIdx = int64_t;

    // This is a very simple data structure for representing graph
    // We only use this structure to sanitize the input file
    // We map the nodes of the graph to integers. Assume the nodes,
    // after sanitizing, lies between 0 and n-1.
struct Graph {
    std::unordered_map<VertexIdx, std::set<VertexIdx>> adjacency_list;
    std::unordered_map<VertexIdx, VertexIdx> degrees; //number of vertices in the graph
};

void add_undirected_edge (Graph g, VertexIdx v1, VertexIdx v2);
VertexIdx count_edges (Graph g);
}
#endif //SUBGRAPHCOUNT_SANITIZE_H
