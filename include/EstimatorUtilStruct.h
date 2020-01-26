//
// Created by Suman Kalyan Bera on 2019-10-09.
//

#ifndef SUBGRAPHCOUNT_ESTIMATORUTILSTRUCT_H
#define SUBGRAPHCOUNT_ESTIMATORUTILSTRUCT_H

#include <string>
#include <vector>
#include <unordered_set>
#include "Graph.h"


using namespace Escape;

struct Parameters {
    std::string filename;
    VertexIdx seed_count = 1;
    std::vector<VertexIdx > seed_vertices;
    EdgeIdx walk_length;
    EdgeIdx subsample_size;
    int no_of_repeat = 1;
    double sparsification_prob = 0.1;
    std::string algo_name = "none!";
    bool print_to_console = false;
    bool print_to_file = true;
    bool edge_count_available = true;
};

// A structure to store an edge with all the relavant information
struct OrderedEdge {
    VertexIdx u;  // the lower degree end point
    VertexIdx v; // the higher degree end point
    EdgeIdx index;   // nbors[index] represents (src, dest)
    VertexIdx degree; // deg(u)
};

// A collection of OrderedEdge -- Used by random walk based algorithms
struct OrderedEdgeCollection {
    EdgeIdx no_of_edges; // No of edges in the collection
    std::vector<OrderedEdge> edge_list; // The list of edges
    std::vector<bool > visited_edge_set; // The set of unique edges in the collection
    std::vector<bool > visited_vertex_set; // The set of unique vertices in the collection
};

// A structure to hold the estimated triangle count and the fraction of the graph seen in the process by an estimator.
struct Estimates {
    double estimate = 0.0;
    double fraction_of_vertices_seen = 1.0;
    double fraction_of_edges_seen = 1.0; // This fraction is with respect to twice the number of edges
};

bool ComparatorByEdgesSeen(Estimates a, Estimates b) {
    return (a.fraction_of_edges_seen < b.fraction_of_edges_seen);
}

bool ComparatorByVerticesSeen(Estimates a, Estimates b) {
    return (a.fraction_of_vertices_seen < b.fraction_of_vertices_seen);
}

// A structure to store the details of vertices and edges observed during the algorithm
struct ObservedGraphStats {
    VertexIdx VerticesSeenInRandomWalk;
    VertexIdx VerticesSeenAsNbors;
    EdgeIdx EdgesSeenInRandomWalk;
    EdgeIdx EdgesSeenAsNbors;
};


// TODO consider merging the below two structure into one
// A structure to store various statistics of an estimator across multiple runs

struct EstimatorStats {
    int no_of_repeats;
    double mean_error_percentage;
    double median_error_percentage;
    double stddev_error_percentage;
    double max_error_percentage;
    double vertices_seen_max_percentage;
    double edges_seen_max_percentage;
};

#endif //SUBGRAPHCOUNT_ESTIMATORUTILSTRUCT_H
