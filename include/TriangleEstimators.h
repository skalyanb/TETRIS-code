//
// Created by Suman Kalyan Bera on 2019-06-17.
//

#ifndef SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H
#define SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H

#include <random>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <string>
#include <ctime>
#include <iterator>


#include <algorithm>
#include "GraphIO.h"
#include "Triadic.h"
#include "Graph.h"
#include "Digraph.h"

using namespace Escape;

/**
 * A structure to keep all the input parameters to various
 triangle estimators.
 * filename: the path to the datafile in escape format
 * seed_count: the number of seed vertices to start random walks from (only for random walk based estimators)
 * walk_length: the length of each random walk (only for random walk based estimators)
 * subsample_size:
 * no_of_repeat: The number of time the estimators are going to be run.
 * sparsification_prob: The sampling probability of the edges in graph sparsification based estimators
 */

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
    std::unordered_set<EdgeIdx> visited_edge_set; // The set of unique edges in the collection
    std::unordered_set<VertexIdx> visited_vertex_set; // The set of unique vertices in the collection
};

// A structure to hold the estimated triangle count and the fraction of the graph seen in the process by an estimator.
struct Estimates {
    double triangle_estimate = 0.0;
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


/**
 * Exactly counting the number of triangles in a graph.
 * This subroutine can be used to count the number of triangles in sparsified graphs
 * @param cg
 * @return
 */
Estimates CountExactTriangles (CGraph *cg)
{
//    printf ("Exactly counting the number of triangles...\n");

    Estimates output;
    cg->sortById();
    CGraph cg_relabel = cg->renameByDegreeOrder();
    cg_relabel.sortById();
    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();

    //printf("Conversion comple. Going to count triangles now. \n");
    TriangleInfo info;
    info = betterWedgeEnumerator(&(dag.outlist));
    output.triangle_estimate = info.total;

    // Free up the memory
    delCGraph(dag.inlist);
    delCGraph(dag.outlist);
    delCGraph(cg_relabel);
    //printf ("  Triangle=%lld ",info.total);
    return output;
}

/**
 * Algorithm:
 * 1. Sample a vertex u.a.r. Call this vertex u.
 * 2. Repeat p many times: // Start a random walk from u.
 * 3.   Sample a neighbor of u u.a.r. Call this vertex v.
 * 4.   Let e= (u,v) and Assume u < v according to degree ordering.
 * 5.   Sample a neighbor w of u u.a.r.
 * 6.   Check for triangle (e,w) where u < v < w according to the degree ordering.
 * 7.   If triangle found, then set Z = d_u, else set Z = 0.
 * 8.   Y += scaling * Z
 * 9.   set u=v (move on to the next step)
 *10. X = Y/p (average over all the steps)
 */

/**
 * There needs to be three functions:
 * 1.   First one does a random walk and collect a long path (collection of edges)
 *      This should be stored in a structure with two elements: edge index and edge degree.
 * 2. Now there would be few variations for comparison.
 *      First variation: all edges collected so far are used for triangle estimation.
 *      Second variation: sample some edges proportional to their degrees and then
 *                        use them for triangle estimation.
 */

OrderedEdgeCollection GetEdgesByRandomWalk(CGraph *cg, Parameters params, std::mt19937 mt) {
    // Get the number of vertices, number of edges, Initialize scaling to number of edges.
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    VertexIdx seed_count = params.seed_count;
    EdgeIdx walk_length = params.walk_length;

    // The list of edges found by the random walk process
    std::vector<OrderedEdge> edge_list;

    // Keep track of the vertices and edges seen so far by the random walk
    std::unordered_set<EdgeIdx> visited_edge_set;
    std::unordered_set<VertexIdx> visited_vertex_set;

    std::vector<EdgeIdx> visited_edge_list;
    std::vector<VertexIdx> visited_vertex_list;


    // Initialize a uniform distribution for generating seed vertices
    //std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);

    // Perform a random walk for seed_count many times
    for (VertexIdx sC = 0; sC < seed_count; sC++) {
        // Pick a random seed vertex
        //VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
        VertexIdx seed = params.seed_vertices[sC];
        VertexIdx parent, child;
        parent = seed;
        // Perform a random walk of length walk_length from the seed vertex
        for (EdgeIdx wL = 0; wL < walk_length; wL++) {
            VertexIdx deg_of_parent = cg->degree(parent); // degree of parent vertex
            // TODO assumption: deg_of_parent is non-zero, add boundary check
            // If the degree of parent vertex is 0, then we attempt seed another vertex
            // We continue this process for certain no of times before giving up and continuing with the run
            int attempt = 0;
            while (deg_of_parent == 0 && attempt < 1000) {
                printf("BAd luck! Isolated vertex! Quitting %d\n", attempt);
                exit(-2);
                //parent = dist_seed_vertex(mt);
                //deg_of_parent = cg->offsets[parent + 1] - cg->offsets[parent];
                //attempt++;
            }

            // Take a step of the random walk
            std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_parent - 1);
            EdgeIdx random_nbor_edge =
                    cg->offsets[parent] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
            child = cg->nbors[random_nbor_edge]; // child is the next vertex on the random walk
            VertexIdx deg_of_child = cg->offsets[child + 1] - cg->offsets[child]; // degree of the child vertex

            // Create the edge info structure
            VertexIdx u, v, low_deg;
            if (deg_of_child < deg_of_parent || (deg_of_child == deg_of_parent && child < parent)) {
                u = child;
                v = parent;
                low_deg = deg_of_child;
            } else {
                u = parent;
                v = child;
                low_deg = deg_of_parent;
            }
            edge_list.push_back(OrderedEdge{u, v, random_nbor_edge, low_deg});

            // Collect the distinct edges and distinct vertices visited so far in a set
//            visited_edge_set.insert(random_nbor_edge);
//            visited_vertex_set.insert(parent);

            visited_edge_list.emplace_back(random_nbor_edge);
            visited_vertex_list.emplace_back(parent);

            // The random walk proceeds with the vertex child
            parent = child;
            // If {u,v} is an isolated edge, then we will stuck here: so find a different seed and continue from there
            if (deg_of_child == 1 && deg_of_parent == 1) {
                printf("BAd luck! Isolated edge. Estimates will be incorrect.\n");
                //parent = dist_seed_vertex(mt);
            }

        }
    }

    // Create return structure for this function.
    // The number of edges produced by the random walk is walk_length;
    std::copy(visited_edge_list.begin(),
              visited_edge_list.end(),
              std::inserter(visited_edge_set, visited_edge_set.end()));

    std::copy(visited_vertex_list.begin(),
              visited_vertex_list.end(),
              std::inserter(visited_vertex_set, visited_vertex_set.end()));

    OrderedEdgeCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
//    printf("Random walk: walk length = %lld ",walk_length);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",visited_edge_set.size(),visited_vertex_set.size());

    return returnEdgeCollection;
}

Estimates SampleByEdgeDegree(CGraph *cg, OrderedEdgeCollection &edge_collection, Parameters params, std::mt19937 mt) {

    // Set up random number generator
    //std::random_device rd;
    //std::mt19937 mt(rd());

    // Collect the degrees of all the edges in a vector so that we can iterate over them
    // The iterator is useful in setting up the discrete_distribution to pick an edge
    // proportionate to its degree
    std::vector<VertexIdx> edge_degree_list(edge_collection.no_of_edges);
    double sum_of_degree = 0.0;
    for (EdgeIdx i = 0; i < edge_collection.no_of_edges; i++) {
        edge_degree_list[i] = edge_collection.edge_list[i].degree;
        sum_of_degree += edge_collection.edge_list[i].degree;
    }

    // std::discrete_distribution is ideal for sampling an edge proportionate to its degree
    // Refer to third answer here: https://stackoverflow.com/questions/1761626/weighted-random-numbers
    std::discrete_distribution<EdgeIdx> wghted_dist_on_edges(edge_degree_list.begin(), edge_degree_list.end());
    //for (double x:wghted_dist_on_edges.probabilities()) std::cout << x << " ";

    // Variables used for estimating the triangle count
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    double X = 0, Y = 0, Z = 0, scaling = edge_collection.no_of_edges;
    Count raw_count = 0; // TODO remove, does not have any purpose

    // Stats about the algorithm: how many new distinct vertices and edges are
    // seen during neighbor sampling process
    ObservedGraphStats obs_graph_stats = {};
    obs_graph_stats.VerticesSeenInRandomWalk = edge_collection.visited_vertex_set.size();
    obs_graph_stats.EdgesSeenInRandomWalk = edge_collection.visited_edge_set.size();

    for (EdgeIdx i = 0; i < params.subsample_size; i++) {

        // Sample an edge proportionate to its degree
        EdgeIdx index = wghted_dist_on_edges(mt);
        OrderedEdge edge = edge_collection.edge_list[index];

        // Sample u.a.r a neighbor w of the edge (of lower degree end point)
        std::uniform_int_distribution<VertexIdx> distNbor(0, edge.degree - 1);
        EdgeIdx random_nbor_edge = cg->offsets[edge.u] + distNbor(mt);
        VertexIdx w = cg->nbors[random_nbor_edge];

        // Update the edges and vertices seen so far data structure
        // TODO to decide whether to pass a pointer or pass a reference for updating the sets
        edge_collection.visited_edge_set.insert(random_nbor_edge);
        edge_collection.visited_vertex_set.insert(w);

        // Now check for a triangle: a triangle is only found if u < v < w and {u,v,w} forms a triangle.
        VertexIdx deg_of_w = cg->offsets[w + 1] - cg->offsets[w];
        VertexIdx deg_of_v = cg->offsets[edge.v + 1] - cg->offsets[edge.v];
        if (cg->isEdgeBinary(w, edge.v) && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && w > edge.v))) {
            EdgeIdx third_edge = cg->getEdgeBinary(edge.v, w);
            edge_collection.visited_edge_set.insert(third_edge);
            //Z = edge.degree;
            Z = 1;
            raw_count++; // Found a triangle
        } else
            Z = 0;
        //Y += scaling * Z;
        Y += Z;
    }
    X = Y / params.subsample_size;
    double sum_degree_estimate = (cg->nEdges / 2.0) * sum_of_degree / edge_collection.no_of_edges;
    X = sum_degree_estimate * X;
    // Create return object and store relevant stats
    Estimates return_estimate = {};
    return_estimate.triangle_estimate = X;
    return_estimate.fraction_of_vertices_seen = edge_collection.visited_vertex_set.size() * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edge_collection.visited_edge_set.size() * 100.0 / m;

    // Some other stats about the algorithm TODO decide what to do with them
    obs_graph_stats.VerticesSeenAsNbors =
            edge_collection.visited_vertex_set.size() - obs_graph_stats.VerticesSeenInRandomWalk;
    obs_graph_stats.EdgesSeenAsNbors = edge_collection.visited_edge_set.size() - obs_graph_stats.EdgesSeenInRandomWalk;

//    printf("**** Sample weighted *****\n");
//    printf("Sum of degree estimate = %lf,actual edges = %lld, err =%lf\n",
//           sum_degree_estimate,cg->nEdges,std::abs(sum_degree_estimate-cg->nEdges)*100/cg->nEdges);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",edge_collection.visited_edge_set.size(),edge_collection.visited_vertex_set.size());
//    printf ("%lf,%lf,%lf\n",X,vertex_fraction,edge_fraction);
    return return_estimate;
}



//EsTRaW
Estimates EstTriByRWandWghtedSampling(CGraph *cg, Parameters params) {
    Estimates output;
    OrderedEdgeCollection randomEdgeCollection;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution
    std::mt19937 mt(rd());

    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    output = SampleByEdgeDegree(cg, randomEdgeCollection, params, mt);
    return output;
}




#endif //SUBGRAPHCOUNT_TRIANGLEESTIMATORS_H
