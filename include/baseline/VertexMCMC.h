//
// Created by Suman Kalyan Bera on 2020-01-22.
//

#ifndef SUBGRAPHCOUNT_VERTEXMCMC_H
#define SUBGRAPHCOUNT_VERTEXMCMC_H


#include <algorithm>

#include "../EstimatorUtilStruct.h"
#include "../TriangleEstimators.h"
#include "../util/RandomWalkUtils.h"


/**
 * Algorithm Details:
 * 1. Start a random walk from a vertex v
 * 2. At each step of the walk, go to a random neighbor of v w.p. min(1,d(v)-1/d(u)-1). Else stay at v.
 * 3. Select two neighbors u.a.r from the current vertex and check if it forms a triangle.
 * 4. Let T be the number of triangles found and W be the number of steps in the random walk.
 * 4. Then transitivity = T/W
 * 4. Triangle Count Estimation = (1/3) * transitivity * \sum_v d(v)* (d(v)-1)/2
 * @param cg
 * @param params
 * @return
 */

Estimates VertexMCMC(CGraph *cg, Parameters params)
{

    /**
     * Setting of various input parameters for executing the algorithm
     */
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // Note that m is double of the number of edges in the graph
    // VertexIdx seed_count = params.seed_count;  // Unused for now. Reserved for future use
    EdgeIdx walk_length = params.walk_length;
    VertexIdx seed = params.seed_vertices[0]; // We are working with only one seed vertex for this project

    /**
    * Set up random number generator
    * Using this random number generator initializize a PRNG: this PRNG is passed along to
    * draw an element from various distributions
    */
    std::random_device rd;
    std::mt19937 mt(rd());

    /**
     * Data structures for book-keeping and local variables
     */
    std::vector <bool > visited_vertex_set(n,false);
    std::vector <bool > visited_edge_set(m,false);
    EdgeIdx count_tri = 0; // This variable counts the number of times we find a triangle
    VertexIdx parent, child,deg_of_parent, deg_of_child;

    /**
     * Algorithm
     */
    parent = seed;
    deg_of_parent = cg->degree(parent); // degree of parent vertex
    if (deg_of_parent == 0)  // If we are stuck with an isolate vertex, we change the seed
        parent = AlternateSeed(cg,mt);

    for (EdgeIdx wL = 0; wL < walk_length; wL++) {  // Perform a random walk of length walk_length from the seed vertex
        deg_of_parent = cg->degree(parent); // degree of parent vertex
        std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_parent - 1);
        EdgeIdx random_nbor_edge = cg->offsets[parent] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
        child = cg->nbors[random_nbor_edge]; // child is the potential next vertex on the random walk
        deg_of_child = cg->degree(child);// degree of the child vertex

        /**
         * Update data structures
         */
        visited_vertex_set[parent] = true;
        visited_vertex_set[child] = true;
        visited_edge_set[random_nbor_edge] = true;

        /**
         * Decide if to take the step or stay put at the same vertex
         */
        bool accept;
        double acceptance_prob;
        if (deg_of_parent == 1)
            accept = true;
        else {
            acceptance_prob = (deg_of_child - 1) * 1.0 / (deg_of_parent - 1);
            std::bernoulli_distribution accpet_dist(acceptance_prob);
            accept = accpet_dist(mt);
        }
        if (accept) {
            parent = child;
            deg_of_parent = deg_of_child;
        }

        /**
         * Sample two nodes from the the neighborhood of parent. Call them u and v. Check if they form a triangle
         */
        std::uniform_int_distribution<VertexIdx> dist_nbor(0, deg_of_parent - 1);
        EdgeIdx random_nbor_edge_u = cg->offsets[parent] + dist_nbor(mt); // TODO: check randomness. same seeding ok?
        EdgeIdx random_nbor_edge_v = cg->offsets[parent] + dist_nbor(mt);
        VertexIdx u = cg->nbors[random_nbor_edge_u];
        VertexIdx v = cg->nbors[random_nbor_edge_v];
        EdgeIdx e = cg->getEdgeBinary(u,v);
        if( u!=v && e != -1) {// parent,u,v forms a triangle
            count_tri++;
        }

        /**
         * Update data structures
         */
        visited_vertex_set[u] = true;
        visited_vertex_set[v] = true;
        visited_edge_set[random_nbor_edge_u] = true;
        visited_edge_set[random_nbor_edge_v] = true;
        if (e!= -1)
            visited_edge_set[e] =true;

    }

    /**
     * Estimate the number of triangles. We need the wedge count for this.
     */
    double wedge_count = 0;
    for (VertexIdx src=0; src < n; src++) {
        VertexIdx deg = cg->degree(src);
        wedge_count += deg * (deg-1) * 0.5;
    }

    double transitivity = 1.0*count_tri / walk_length;
    double triangleEstimate = transitivity * wedge_count / 3.0 ;

    /**
     * Populate the return structure
     */
    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;

    VertexIdx edges_seen = std::count(visited_edge_set.begin(),visited_edge_set.end(),true);
    VertexIdx vertices_seen = std::count(visited_vertex_set.begin(),visited_vertex_set.end(), true);

    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / n;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / m;

    return return_estimate;
}
#endif //SUBGRAPHCOUNT_VERTEXMCMC_H
