//
// Created by Suman Kalyan Bera on 2019-07-19.
//

// Takes as input a file with edges list as pairs of strings, and creates
// the proper input for our library.
// Usage: sanitize.out out_directory input_filename
// where input_filename is expected to be a relative path to the input file
// and out_directory is where the .edges file will be produced in proper format.
// This format has the first line with the number of nodes and edges.
// Each line has a distinct edge with node labels as ints starting from 0.

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>
#include <unordered_map>

#include "sanitize.h"

using namespace sanitize;

void sanitize::add_undirected_edge(sanitize::Graph g, sanitize::VertexIdx v1, sanitize::VertexIdx v2)
{
    // If self loop, then do nothing
    if (v1 == v2) {
    return;
    }

    // Find if v1 is already present in the graph
    auto it1 = g.adjacency_list.find (v1);

    if (it1 != g.adjacency_list.end()) {        // Check if v1 is part of the graph
        // Retrieve the adjacency list of v1. Since the adjacency list is maintained as a set,
        // we can safely insert v2 in this without worrying about duplicate entries.
        auto ret = it1->second.emplace(v2);
        // If v2 is not already present, then emplace returns true as the second element is the pair
        if (ret.second) {
            g.degrees[v1] += 1;
        }
    }
    else {                                      // v1 is not part of the graph yet
        std::set<VertexIdx> nbr = {v2};         // Initialize the neighbor list for v1
        g.adjacency_list[v1] = nbr;             // Insert v1` in the graph
        g.degrees[v1] = 1;                      // Update the degree of the graph
    }

    // Now repeat the above process for the vertex v2
    auto it2 = g.adjacency_list.find (v2);

    if (it2 != g.adjacency_list.end()) {        // Check if v2 is part of the graph
        auto ret = it2->second.emplace(v1);
        if (ret.second) {
            g.degrees[v2] += 1;
        }
    }
    else {                                      // v2 is not part of the graph yet
        std::set<VertexIdx > nbr = {v1};
        g.adjacency_list[v2] = nbr;
        g.degrees[v2] = 1;
    }
}

VertexIdx sanitize::count_edges (Graph g)
{
    VertexIdx m = 0;
    for (auto it=g.degrees.begin(); it!= g.degrees.end(); it++)
        m += it->second;
    return m/2;
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cout << "Usage: sanitize.out out_directory input_filename" << std::endl;
        return 0;
    }
    // Open the input file in read only mode
    std::ifstream input_file(argv[2]);
    if (!input_file.is_open()) {
        std::cout << "Error in opening input file." << std::endl;
        return 0;
    }

    // Create a map that will map every vertex id the in the input file to an integer in the output file
    //std::map <sanitize::VertexIdx, sanitize::VertexIdx > dict;
    std::unordered_map <VertexIdx, VertexIdx> dict;

    std::string line;
    VertexIdx index = 0;

    // Create an empty graph
    Graph g;

    while (std::getline (input_file,line)) {

        // Note: here the assumption is that the delimiter is whitespace
        // If we need to deal with some other delimiter, then we need to
        // incorporate a separate split function
        std::istringstream iss (line);
        std::vector <std::string> token = {std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>{}};

        // If the first character is # or % then ignore that line
        // If there are not exactly two entries, then also ignore the line
        if (token[0] == "#" || token[0] == "%" || token.size()!=2)
            continue;

        // Remove self loops
        if ( token[0] == token[1])
            continue;

        // We are assuming the input file has integers as vertex id
        std::string::size_type sz = 0;   // alias of size_t
        VertexIdx node1 = std::stoll(token[0],&sz,0);
        VertexIdx node2 = std::stoll(token[1],&sz,0);

        // Find if the vertices are already present in the dictionary
        auto it1 = dict.find (node1);
        auto it2 = dict.find (node1);
        // If the first vertex is present then replace the node1 with the corresponding node id;
        // Otherwise create a new node for this vertex.
        if (it1!= dict.end())
            node1 = it1->second;
        else {
            dict[node1] = index;
            node1 = index;
            index++;
        }
        // Repeat the same for the second node
        if (it2!= dict.end())
            node2 = it2->second;
        else {
            dict[node2] = index;
            node2 = index;
            index++;
        }
        // Add the edge (node1,node2) to the graph
        // Note that during this process, we did not check if this is an duplicate edge
        // Duplicate edges will be discarded during the execution of the following function
        add_undirected_edge(g,node1,node2);
    }
    input_file.close();
    VertexIdx n = g.degrees.size(); // number of vertices
    VertexIdx m = count_edges(g);   // number of edges

    // Now we create an output file with the same name as input file in the output directory


    // Finally we write the graph: there is no duplicate edges in the output file, every edge appears exactly once

//    for
//    for node in G.vertices:
//    for nbr in G.adj_list[node]:
//    if node == n:
//    print node, "node out of range"
//    if nbr == n:
//    print nbr, "nbr out of range"
//    if node < nbr:
//    f_output.write(str(node)+' '+str(nbr)+'\n')
}