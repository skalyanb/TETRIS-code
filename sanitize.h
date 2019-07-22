//
// Created by Suman Kalyan Bera on 2019-07-19.
//

#ifndef SUBGRAPHCOUNT_SANITIZE_H
#define SUBGRAPHCOUNT_SANITIZE_H

#include <cstdlib>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <vector>

namespace sanitize
{
    using VertexIdx = int32_t;

    // This is a very simple data structure for representing graph
    // We only use this structure to sanitize the input file
    // We map the nodes of the graph to integers. Assume the nodes,
    // after sanitizing, lies between 0 and n-1.
    //bool comp (VertexIdx lhs, VertexIdx rhs) {return lhs<rhs;}

    struct Graph {
    std::map<VertexIdx, std::set<VertexIdx>> adjacency_list;
    //std::map<std::string, std::set<VertexIdx>> adjacency_list;
    //std::unordered_map<VertexIdx, std::unordered_set<VertexIdx>> adjacency_list;
    //std::unordered_map<VertexIdx, std::set<VertexIdx>> adjacency_list;
    //std::map<VertexIdx, VertexIdx> degrees; //number of vertices in the graph
    std::unordered_map<VertexIdx, VertexIdx> degrees; //number of vertices in the graph
    };

    void add_undirected_edge (Graph &g, VertexIdx v1, VertexIdx v2) {
        // If self loop, then do nothing
        if (v1 == v2) {
            return;
        }

        // Find if v1 is already present in the graph
        auto it1 = g.adjacency_list.find(v1);

        if (it1 != g.adjacency_list.end()) {        // Check if v1 is part of the graph
            // Retrieve the adjacency list of v1. Since the adjacency list is maintained as a set,
            // we can safely insert v2 in this without worrying about duplicate entries.
            auto ret = it1->second.emplace(v2);
            // If v2 is not already present, then emplace returns true as the second element is the pair
            if (ret.second) {
                g.degrees[v1] += 1;
            }
        } else {                                      // v1 is not part of the graph yet
            //std::unordered_set<VertexIdx> nbr = {v2};         // Initialize the neighbor list for v1
            std::set<VertexIdx> nbr = {v2};         // Initialize the neighbor list for v1
            g.adjacency_list[v1] = nbr;             // Insert v1` in the graph
            g.degrees[v1] = 1;                      // Update the degree of the graph
        }

        // Now repeat the above process for the vertex v2
        auto it2 = g.adjacency_list.find(v2);

        if (it2 != g.adjacency_list.end()) {        // Check if v2 is part of the graph
            auto ret = it2->second.emplace(v1);
            if (ret.second) {
                g.degrees[v2] += 1;
            }
        } else {                                      // v2 is not part of the graph yet
            //std::unordered_set<VertexIdx > nbr = {v1};
            std::set<VertexIdx> nbr = {v1};
            g.adjacency_list[v2] = nbr;
            g.degrees[v2] = 1;
        }
    }
    VertexIdx count_edges (Graph g)
    {
        VertexIdx m = 0;
        for (auto it=g.degrees.begin(); it!= g.degrees.end(); it++)
            m += it->second;
        return m/2;
    }
}
#endif //SUBGRAPHCOUNT_SANITIZE_H
