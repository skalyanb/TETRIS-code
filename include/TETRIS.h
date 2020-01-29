//
// Created by Suman Kalyan Bera on 2020-01-25.
//

#ifndef SUBGRAPHCOUNT_TETRIS_H
#define SUBGRAPHCOUNT_TETRIS_H

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
#include "EstimatorUtilStruct.h"
#include "EstimateEdgeCount.h"
#include "util/RandomWalkUtils.h"


using namespace Escape;

OrderedEdgeCollection GetEdgesByRandomWalk(CGraph *cg, Parameters params, std::mt19937 mt) {
    /**
     * Setting of various input parameters for executing the algorithm
     */
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges;
    VertexIdx seed_count = params.seed_count; // We are working with only one seed vertex for this project
    EdgeIdx walk_length = params.walk_length;

    /**
     * Data structures for book-keeping and storing the edges found by the random walk
     */
    std::vector<OrderedEdge> edge_list;
    std::vector<bool > visited_edge_set(m,false);
    std::vector<bool > visited_vertex_set(n,false);
    VertexIdx no_of_query = 0;  // Counts the number of query made by the algorithm

    /**
     * Local variables
     */
    VertexIdx parent, child,deg_of_parent, deg_of_child;

    for (VertexIdx sC = 0; sC < seed_count; sC++) {

         VertexIdx seed = params.seed_vertices[sC]; // Random seed vertex is already in the params structure.
            parent = seed;
            deg_of_parent = cg->degree(parent); // degree of parent vertex
            if (deg_of_parent == 0)  // If we are stuck with an isolate vertex, we change the seed
                parent = AlternateSeed(cg,mt);

            for (EdgeIdx wL = 0; wL < walk_length; wL++) { // Perform a random walk of length walk_length from the seed vertex and collect the edges
                deg_of_parent = cg->degree(parent); // degree of parent vertex
                std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_parent - 1);
                EdgeIdx random_nbor_edge = cg->offsets[parent] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
                no_of_query++; // One query to the uniform random neighbor oracle
                child = cg->nbors[random_nbor_edge]; // child is the next vertex on the random walk
                deg_of_child = cg->degree(child); //degree of child vertex

                /**
                 * Update the edge collection data structure with relevant information about the edge
                 */
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

                /**
                 * Update the book-keeping data structure with the visited edge and vertices
                 */
                visited_edge_set[random_nbor_edge] = true;
                visited_vertex_set[parent] = true;
                parent = child; // The random walk proceeds with the vertex child
            }
        }

    OrderedEdgeCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set,no_of_query};
    return returnEdgeCollection;
}

/**
 * In this function, we simulate degree based sampling from the set of collected edges and estimate the triangle count based on that.
 * @param cg
 * @param edge_collection
 * @param params
 * @param mt
 * @param edge_estimate
 * @return
 */
Estimates SampleByEdgeDegree(CGraph *cg, OrderedEdgeCollection &edge_collection, Parameters params, std::mt19937 mt, double edge_estimate) {

    /**
     * Collect the degrees of all the edges in a vector so that we can iterate over them
     * The iterator is useful in setting up the discrete_distribution to pick an edge
     * proportionate to its degree
     */
    std::vector<VertexIdx> edge_degree_list(edge_collection.no_of_edges);
    double sum_of_degree = 0.0;
    for (EdgeIdx i = 0; i < edge_collection.no_of_edges; i++) {
        edge_degree_list[i] = edge_collection.edge_list[i].degree;
        sum_of_degree += edge_collection.edge_list[i].degree;
    }

    /**
     * We next set up a distribution for degree-based sampling.
     * std::discrete_distribution is ideal for sampling an edge proportionate to its degree
     * Refer to third answer here: https://stackoverflow.com/questions/1761626/weighted-random-numbers
     */
    std::discrete_distribution<EdgeIdx> wghted_dist_on_edges(edge_degree_list.begin(), edge_degree_list.end());

    /**
     * Local variables -- estimators
     */
    double X = 0, Y = 0, Z = 0;

    for (EdgeIdx i = 0; i < params.subsample_size; i++) {
        EdgeIdx index = wghted_dist_on_edges(mt); // Sample an edge proportionate to its degree
        OrderedEdge edge = edge_collection.edge_list[index]; // Retrieve the edge information

        /**
         * Sample u.a.r a neighbor w of the edge (of lower degree end point)
         */
        std::uniform_int_distribution<VertexIdx> distNbor(0, edge.degree - 1);
        EdgeIdx random_nbor_edge = cg->offsets[edge.u] + distNbor(mt);
        edge_collection.no_of_query++; // One query to the uniform random neighbor oracle
        VertexIdx w = cg->nbors[random_nbor_edge];
        /**
         * Update the edges and vertices seen so far data structure
         */
        edge_collection.visited_edge_set[random_nbor_edge] = true;
        edge_collection.visited_vertex_set[w] = true;
        /**
         * Now check for a triangle: a triangle is only found if u < v < w and {u,v,w} forms a triangle.
         */
        VertexIdx deg_of_w = cg->degree(w);
        VertexIdx deg_of_v = cg->degree(edge.v);
        edge_collection.no_of_query++; // One query to the edge query in the if condition
        if (cg->isEdgeBinary(w, edge.v) && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && w > edge.v))) {
            EdgeIdx third_edge = cg->getEdgeBinary(edge.v, w);
            edge_collection.visited_edge_set[third_edge] = true;
            Z = 1;
        } else
            Z = 0;
        Y += Z;
    }
    X = Y / params.subsample_size;
    double sum_degree_estimate = (edge_estimate / 2.0) * sum_of_degree / edge_collection.no_of_edges;
    X = sum_degree_estimate * X;

    /**
     * Create return object and store relevant stats
     */
    Estimates return_estimate = {};
    return_estimate.estimate = X;
    VertexIdx vertices_seen = std::count(edge_collection.visited_vertex_set.begin(),edge_collection.visited_vertex_set.end(),true);
    EdgeIdx edges_seen = std::count(edge_collection.visited_edge_set.begin(),edge_collection.visited_edge_set.end(),true);
    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / cg->nEdges;
    return_estimate.query_complexity = edge_collection.no_of_query * 100.0 / cg->nEdges;

    return return_estimate;
}

Estimates TETRIS(CGraph *cg, Parameters params) {

    /**
     * Set up local variables. skip parameter is used in estimating the number of edges in the graph.
     */
    Estimates output, edge_count_output;
    OrderedEdgeCollection randomEdgeCollection;
    int skip = 25;

    /**
    * Set up random number generator
    * Using this random number generator initializize a PRNG: this PRNG is passed along to
    * draw an element from various distributions
    */
    std::random_device rd;
    std::mt19937 mt(rd());

    /**
     * Perform a random walk through the graph, and collect the edges in the process
     */
    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);

    /**
     * Depending on whether the total number of edges in the graphs is available or not, compute it.
     */
    if (!params.normalization_count_available)
        edge_count_output = EstimateEdgeCount (cg,randomEdgeCollection, params, skip);
    else
        edge_count_output.estimate = cg->nEdges *1.0;

    /**
     * Perform a weighted sampling and estimate the triangle count
     */
    output = SampleByEdgeDegree(cg, randomEdgeCollection, params, mt, edge_count_output.estimate);
    return output;
}


#endif //SUBGRAPHCOUNT_TETRIS_H
