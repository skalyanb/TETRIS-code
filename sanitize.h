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
    using VertexIdx = int64_t;
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

}
#endif //SUBGRAPHCOUNT_SANITIZE_H
