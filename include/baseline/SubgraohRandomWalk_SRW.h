//
// Created by Suman Kalyan Bera on 2020-01-24.
//

#ifndef SUBGRAPHCOUNT_SUBGRAOHRANDOMWALK_SRW_H
#define SUBGRAPHCOUNT_SUBGRAOHRANDOMWALK_SRW_H

#include <random>
#include <vector>

#include "../Graph.h"
#include "../EstimatorUtilStruct.h"
#include "../util/RandomWalkUtils.h"
#include "../EstimateEdgeCount.h"

/**
 * Algorithm Details:
 * 1. Start a random walk from a vertex v1
 * 2. Take 3 steps and assume X0 = (v1,v2,v3) be the vertices visited in that order.
 * 3. Perform a random walk of given length.
 * 4. At each step of the walk, consider the last 3 vertices visited by the walk.
 * 5a. If it forms a triangle, and CSS flag is not set, then add m*d(v2)/3.
 * 5b. If it forms a triangle, and CSS flag is set, then add m/(1/d(v1)+1/(dv2)+1/d(v3)).
 * 6. Finally, take the average over the length of the random walk.
 * Note that the factor m/3 can be multiplied at the end of the algorithm
 * The parameter NB decides if the random walk is non-backtracking.
 * If it is true, then v1-v2-v1 becomes an invalid state.
 * We sample another random neighbor in that case.
 * Also, in final estimation, d'(u)=d(u)-1 is used.
 * We are assuming NB is set to true only when CSS is also set to true
 * @param cg
 * @param params
 * @return
 */

Estimates SRW1(CGraph *cg, Parameters params)
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
     * Data structures for book-keeping and local variables.
     * We do not really need the OrderedEdge vector for SRW.
     * However, if m is not available, and needs to be estimated based on the random walk,
     * then we need the collection of edges.
     */
    std::vector<OrderedEdge> edge_list;
    std::vector <bool > visited_vertex_set(n,false);
    std::vector <bool > visited_edge_set(m,false);
    VertexIdx no_of_query = 0;  // Counts the number of query made by the algorithm
    double X = 0; // This variable keeps a running count of the estimator
    VertexIdx v1, v2, v3, deg_of_v1, deg_of_v2, deg_of_v3;
    EdgeIdx random_nbor_edge;

    /**
     * Algorithm Initialization: Take two steps of the random walk and collect v1 and v2.
     */

    v1 = seed;
    deg_of_v1 = cg->degree(v1); // degree of parent vertex
    if (deg_of_v1 == 0)  // If we are stuck with an isolate vertex, we change the seed
        v1 = AlternateSeed(cg,mt);
    deg_of_v1 = cg->degree(v1);
    std::uniform_int_distribution<VertexIdx> dist_parent(0, deg_of_v1 - 1);
    no_of_query++; // One query to the uniform random neighbor oracle
    random_nbor_edge = cg->offsets[v1] + dist_parent(mt); // TODO: check randomness. same seeding ok?
    v2 = cg->nbors[random_nbor_edge]; // child is the potential next vertex on the random walk
    deg_of_v2 = cg->degree(v2);// degree of the child vertex

    /**
    * Update data structures
    */
    edge_list.push_back(OrderedEdge{v1, v2, random_nbor_edge, deg_of_v1}); // Note that we do not care of the ordering of the edge.
    visited_vertex_set[v1] = true;
    visited_vertex_set[v2] = true;
    visited_edge_set[random_nbor_edge] = true;


    for (EdgeIdx wL = 0; wL < walk_length; wL++) {  // Perform a random walk of length walk_length from the vertex v2
        std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_v2 - 1); // v2 is the current last visited vertex on the random walk
        random_nbor_edge = cg->offsets[v2] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
        no_of_query++; // One query to the uniform random neighbor oracle
        v3 = cg->nbors[random_nbor_edge]; // v3 is the  new vertex on the random walk
        deg_of_v3 = cg->degree(v3);// degree of the vertex v3

        /**
         * The parameter NB decides if the random walk is non-backtracking.
         * If it is true, then v1-v2-v1 becomes an invalid state.
         * We sample another random neighbor in that case.
         * Also, in final estimation, d'(u)=d(u)-1 is used.
         */
        if (params.NB) {
            while (v1 == v3) { // We reached an invalid state, resample.
                random_nbor_edge = cg->offsets[v2] + dist_parent_nbor(mt); // TODO: check randomness. same seeding ok?
                no_of_query++; // One query to the uniform random neighbor oracle
                v3 = cg->nbors[random_nbor_edge]; // v3 is the  new vertex on the random walk
                deg_of_v3 = cg->degree(v3);// degree of the vertex v3
            }
        }

        /**
         * Update data structures
         */
        edge_list.push_back(OrderedEdge{v2, v3, random_nbor_edge, deg_of_v2}); // Note that we do not care of the ordering of the edge.
        visited_vertex_set[v3] = true;
        visited_edge_set[random_nbor_edge] = true;

        /**
         * Check if v1, v2 and v3 forms a triangle
         * Note that (v1,v2) and (v2,v3) edges are present.
         * So we check if (v1,v3) edge is present and if all v1,v2 and v3 are distinct
         */
        EdgeIdx e = cg->getEdgeBinary(v1,v3);
        no_of_query++; // One query to the edge query oracle
        if (e!= -1 && v1!=v2 && v2!=v3 && v1!=v3) {
            if (params.CSS && params.NB) {  // Both CSS and NB are true
                double denom = (1.0 / (deg_of_v1-1)) +
                                (1.0 /(deg_of_v2-1)) +
                                (1.0 /(deg_of_v3-1));
                X += 1.0 / denom;
            }
            else if (params.CSS){ // Only CSS is true
                double denom = (1.0 / deg_of_v1) + (1.0 /deg_of_v2) + (1.0 /deg_of_v3);
                X += 1.0 / denom;
            }
            else if (params.NB) // Only NB is true
                X+= deg_of_v2 - 1;
            else                // Both CSS and NB are false
                X += deg_of_v2;


        }

        /**
         * Update data structures
         */
        if (e!= -1)
            visited_edge_set[e] =true;
        /**
         * slide the window of history!
         */
        v1 = v2;
        deg_of_v1 = deg_of_v2;
        v2 = v3;
        deg_of_v2 = deg_of_v3;
    }

    double edge_estimate = 0.0;
    /**
     * Depending on whether the total number of edges in the graphs is available or not, compute it.
     */
    if (!params.normalization_count_available) {
        int skip = 100;
        OrderedEdgeCollection randomEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
        Estimates edge_count_output = EstimateEdgeCount(cg, randomEdgeCollection, params, skip);
        edge_estimate = edge_count_output.estimate;
    }
    else
        edge_estimate = cg->nEdges *1.0;

    /**
     * Estimate the number of triangles.
     */
    double triangleEstimate = 0.0;
    if (params.CSS)
        triangleEstimate = (1.0/walk_length) * X * (edge_estimate / 2.0) ; // Note that m is actually twice the number of edges
    else
        triangleEstimate = (1.0/walk_length) * X * (edge_estimate/2.0) * (1/ 3.0) ; // Note that m is actually twice the number of edges

    /**
     * Populate the return structure
     */
    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;

    VertexIdx edges_seen = std::count(visited_edge_set.begin(),visited_edge_set.end(),true);
    VertexIdx vertices_seen = std::count(visited_vertex_set.begin(),visited_vertex_set.end(), true);

    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / n;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / m;
    return_estimate.query_complexity = no_of_query * 100.0 /m;

    return return_estimate;
}


#endif //SUBGRAPHCOUNT_SUBGRAOHRANDOMWALK_SRW_H
