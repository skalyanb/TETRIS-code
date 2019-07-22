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
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility> // pair

#include "sanitize.h"

using namespace sanitize;

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1, T2> const &pair) const
    {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);

        return h1 ^ h2;
    }
};

bool comp (std::pair<VertexIdx ,VertexIdx > first, std::pair<VertexIdx ,VertexIdx > second)
{
    if (first.first == second.first )
        return first.second < second.second;
    else
        return first.first < second.first;
}

bool equal (std::pair<VertexIdx ,VertexIdx > first, std::pair<VertexIdx ,VertexIdx > second)
{
    return (first.first == second.first && first.second == second.second );
}

int sanitize_direct (char* argv[])
{
    // Open the input file in read only mode
    std::ifstream input_file(argv[2]);
    if (!input_file.is_open()) {
        std::cout << "Error in opening input file." << std::endl;
        return 0;
    }
    // Create a map that will map every vertex id the in the input file to an integer in the output file
    std::unordered_map <VertexIdx, VertexIdx> dict;
    dict.reserve(100000000); // 100M

    // Create an unordered set for fast storing of edges as a pair of vertices
    std::unordered_set<std::pair<VertexIdx ,VertexIdx >, pair_hash> edgeList;
    //edgeList.reserve(500000000); // 500M

    // Create a vector for storing pair of vertices, later we will sort and remove duplicates from this list
    std::vector<std::pair<VertexIdx ,VertexIdx >> edgeListDup;
    //edgeList.reserve(500000000); // 500M

    std::string line;
    VertexIdx index = 0;

    // Debug help messages and parameters
    long long int line_no = 0; // just for progress checking pupose -- remove
    // Check the capacity of the hash table
    std::cout << "max_size = " << dict.max_size() << ";";
//              << " max_bucket_count = " << dict.max_bucket_count() << ";"
//              << "max_load_factor = " << dict.max_load_factor() << std::endl;
    // Check the capacity of the edge list etc
//    std::cout << "max_size = " << edgeList.max_size() << ";"
//              << " max_bucket_count = " << edgeList.max_bucket_count() << ";"
//              << "max_load_factor = " << edgeList.max_load_factor() << std::endl;

    while (std::getline (input_file,line)) {
        // Note: here the assumption is that the delimiter is whitespace
        // If we need to deal with some other delimiter, then we need to
        // incorporate a separate split function
        std::istringstream iss (line);
        std::vector <std::string> token = {std::istream_iterator<std::string>{iss},
                                           std::istream_iterator<std::string>{}};

        // If the first character is # or % then ignore that line
        // If there are not exactly two entries, then also ignore the line
        if (token[0][0] == '#' || token[0][0] == '%' || token.size()!=2)
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
        auto it2 = dict.find (node2);
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
        //Add the edge (node1,node2) as a pair to the edgeList. Ensure node1 < node2
        if (node1 < node2)
//            edgeList.emplace(std::make_pair(node1,node2));
            edgeListDup.emplace_back(std::make_pair(node1,node2));
        else
//            edgeList.emplace(std::make_pair(node2,node1));
            edgeListDup.emplace_back(std::make_pair(node2,node1));

        // Debug help messages and parameters
        line_no ++;
        if (line_no % 1000000 == 0){
            std::cout << "At line no " << line_no << std::endl;
            std::cout << "current size = " << dict.size() << ";";
//                      << " current_bucket_count = " << dict.bucket_count() << ";"
//                      << "current_load_factor = " << dict.load_factor() << std::endl;
            // Check the capacity of the edge list etc
//            std::cout << "current _size = " << edgeList.size() << ";"
//                      << " current_bucket_count = " << edgeList.bucket_count() << ";"
//                      << "current_load_factor = " << edgeList.load_factor() << std::endl;
        }
    }
    input_file.close();
    std::sort(edgeListDup.begin(),edgeListDup.end(),comp);
    edgeListDup.erase(std::unique(edgeListDup.begin(),edgeListDup.end(),equal),edgeListDup.end());

    std::cout << "File reading complete." << std::endl;
    // Number of vertices and edge
    VertexIdx n = dict.size();
    VertexIdx m = edgeListDup.size();

    // Sort the edgeList according to vertexIdx
    //std::vector<std::pair<VertexIdx ,VertexIdx >> edgeVector(edgeList.begin(),edgeList.end());
    //std::sort (edgeVector.begin(), edgeVector.end(), comp);

    // Now we create an output file with the same name as input file in the output directory
    std::string filename_with_path(argv[2]);
    std::string outfile_dir(argv[1]);
    size_t ind_begin = filename_with_path.find_last_of("/\\");
    size_t ind_end = filename_with_path.find_last_of(".");
    std::string out_file_path = outfile_dir + filename_with_path.substr(ind_begin + 1, ind_end - ind_begin - 1)
                                + ".edges";
    std::ofstream output_file(out_file_path);
    if (!output_file.is_open()) {
        std::cout << "File writing error." << std::endl;
        return 0;
    }

    // First line gets the number of vertices and number of edges: n m
    output_file << n << " " << m << "\n";
    for (auto element : edgeListDup)
        output_file << element.first << " " << element.second << "\n";

    // Finally we write the graph: there is no duplicate edges in the output file, every edge appears exactly once

    output_file.close();
    return 0;
}


int sanitize_with_graph (char *argv[])
{
    // Open the input file in read only mode
    std::ifstream input_file(argv[2]);
    if (!input_file.is_open()) {
        std::cout << "Error in opening input file." << std::endl;
        return 0;
    }

    // Create a map that will map every vertex id the in the input file to an integer in the output file
    //std::map <VertexIdx, VertexIdx > dict;
    std::unordered_map <VertexIdx, VertexIdx> dict;
    dict.reserve(500000000); // 500M

    std::string line;
    VertexIdx index = 0;

    // Create an empty graph
    Graph g;
    //g.degrees.reserve(500000000); // 500M

    // Debug help messages and parameters
    long long int line_no = 0; // just for progress checking pupose -- remove
    // Check the capacity of the hash table
    std::cout << "max_size = " << dict.max_size() << ";"
              << " max_bucket_count = " << dict.max_bucket_count() << ";"
              << "max_load_factor = " << dict.max_load_factor() << std::endl;
    while (std::getline (input_file,line)) {
        // Note: here the assumption is that the delimiter is whitespace
        // If we need to deal with some other delimiter, then we need to
        // incorporate a separate split function
        std::istringstream iss (line);
        std::vector <std::string> token = {std::istream_iterator<std::string>{iss},
                                           std::istream_iterator<std::string>{}};

        // If the first character is # or % then ignore that line
        // If there are not exactly two entries, then also ignore the line
        if (token[0][0] == '#' || token[0][0] == '%' || token.size()!=2)
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
        auto it2 = dict.find (node2);
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

        line_no ++;
        if (line_no % 10000000 == 0){
            std::cout << "At line no " << line_no << std::endl;
            std::cout << "current size = " << dict.size() << ";";
            std::cout << " current_bucket_count = " << dict.bucket_count() << ";";
            std::cout << "current_load_factor = " << dict.load_factor() << std::endl;
            std::cout << "current size of adjacency list = " << g.adjacency_list.size() << std::endl;
        }
    }
    input_file.close();
    VertexIdx n = g.degrees.size();
    VertexIdx m = count_edges(g);

    // Now we create an output file with the same name as input file in the output directory
    std::string filename_with_path(argv[2]);
    std::string outfile_dir(argv[1]);
    size_t ind_begin = filename_with_path.find_last_of("/\\");
    size_t ind_end = filename_with_path.find_last_of(".");
    std::string out_file_path = outfile_dir + filename_with_path.substr(ind_begin + 1, ind_end - ind_begin - 1)
                                + ".edges";
    std::ofstream output_file(out_file_path);
    if (!output_file.is_open()) {
        std::cout << "File writing error." << std::endl;
        return 0;
    }

    // First line gets the number of vertices and number of edges: n m
    output_file << n << " " << m << "\n";

    // Finally we write the graph: there is no duplicate edges in the output file, every edge appears exactly once
    for (auto node_it = g.adjacency_list.begin(); node_it!= g.adjacency_list.end(); node_it ++) {
        for (auto nbr_it = node_it->second.begin(); nbr_it!= node_it->second.end(); nbr_it++){
            if (node_it->first < *nbr_it ) {
                output_file << node_it->first << " " << *nbr_it << "\n";
            }
        }
    }
    output_file.close();
    return 0;
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cout << "Usage: sanitize.out out_directory input_filename" << std::endl;
        return 0;
    }
    return sanitize_direct(argv);
    //return sanitize_with_graph(argv);
}