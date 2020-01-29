//
// Created by Suman Kalyan Bera on 2020-01-23.
//

#ifndef SUBGRAPHCOUNT_RANDOMWALKUTILS_H
#define SUBGRAPHCOUNT_RANDOMWALKUTILS_H

#include <random>

#include "../Graph.h"

VertexIdx AlternateSeed (CGraph* cg, std::mt19937 mt)
{
    int attempt = 0;
    VertexIdx deg_of_parent = 0, parent = 0;
    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, cg->nVertices - 1); // Initialize a uniform distribution for generating seed vertices
    while (deg_of_parent == 0 && attempt < 1000) {
        printf("BAd luck! Trying again %d\n", attempt);
        parent = dist_seed_vertex(mt);
        deg_of_parent = cg->degree(parent);
        attempt++;
    }
    return parent;
}


OrderedEdgeCollection GetEdgesByRandomWalk_2(CGraph *cg, Parameters params, std::mt19937 mt) {
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
//            VertexIdx u, v, low_deg;
//            if (deg_of_child < deg_of_parent || (deg_of_child == deg_of_parent && child < parent)) {
//                u = child;
//                v = parent;
//                low_deg = deg_of_child;
//            } else {
//                u = parent;
//                v = child;
//                low_deg = deg_of_parent;
//            }
            edge_list.push_back(OrderedEdge{parent, child, random_nbor_edge, deg_of_parent});

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

#endif //SUBGRAPHCOUNT_RANDOMWALKUTILS_H
