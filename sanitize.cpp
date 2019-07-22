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
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility> // pair

#include "sanitize.h"

using namespace sanitize;

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

    // Create a vector for storing pair of vertices, later we will sort and remove duplicates from this list
    // For fast processing, we will simply read an edge and put it at the end of edgeList.
    std::vector<std::pair<VertexIdx ,VertexIdx >> edgeList;

    std::string line;
    VertexIdx index = 0;

    // Debug help messages and parameters
    long long int line_no = 0; // just for progress checking progress of code -- remove
    // Check the capacity of the hash table
//    std::cout << "max_size = " << dict.max_size() << ";";
//              << " max_bucket_count = " << dict.max_bucket_count() << ";"
//              << "max_load_factor = " << dict.max_load_factor() << std::endl;

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
            edgeList.emplace_back(std::make_pair(node1,node2));
        else
            edgeList.emplace_back(std::make_pair(node2,node1));

        // Debug help messages and parameters
        line_no ++;
        if (line_no % 1000000 == 0){
            std::cout << "At line no " << line_no << std::endl;
//            std::cout << "current size = " << dict.size() << ";";
        }
    }
    input_file.close();
    std::cout << "File reading complete." << std::endl;
    // Sort the edgeList and then remove the duplicate by calling unique function. Finally erase the duplicate.
    std::sort(edgeList.begin(),edgeList.end(),comp);
    edgeList.erase(std::unique(edgeList.begin(),edgeList.end(),equal),edgeList.end());

    // Number of vertices and edge
    VertexIdx n = dict.size();
    VertexIdx m = edgeList.size();

    // Now we create an output file with the same name as input file in the output directory
    std::string filename_with_path(argv[2]);
    std::string outfile_dir(argv[1]);
    size_t ind_begin = filename_with_path.find_last_of("/\\");
    size_t ind_end = filename_with_path.find_last_of('.');
    std::string out_file_path = outfile_dir + filename_with_path.substr(ind_begin + 1, ind_end - ind_begin - 1)
                                + ".edges";
    std::ofstream output_file(out_file_path);
    if (!output_file.is_open()) {
        std::cout << "File writing error." << std::endl;
        return 0;
    }

    // First line gets the number of vertices and number of edges: n m
    output_file << n << " " << m << "\n";
    for (auto element : edgeList)
        output_file << element.first << " " << element.second << "\n";

    output_file.close();
    return 0;
}



int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cout << "Usage: sanitize.out out_directory input_filename" << std::endl;
        return 0;
    }
    return sanitize_direct(argv);
}