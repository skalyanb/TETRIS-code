//
// Created by Suman Kalyan Bera on 2020-01-28.
//

#ifndef SUBGRAPHCOUNT_SERWC_H
#define SUBGRAPHCOUNT_SERWC_H

#include <set>
#include <unordered_map>

#include "../EstimatorUtilStruct.h"
#include "../TriangleEstimators.h"
#include "../util/RandomWalkUtils.h"
#include "../util/TrianglleCountUtil.h"

/**
 * Algorithm Details:
 * 1. Start a random walk from a vertex v
 * 2. Perform a random walk of given length.
 * 3. At each step of the walk, count the number of triangles incident on the edge.
 * 4. Finally, take the average over the length of the random walk.
 * 5. Scale up by m/3 where m is the number of edges.
 * @param cg
 * @param params
 * @return
 */

Estimates SERWC(CGraph *cg, Parameters params)
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

    VertexIdx parent, child,deg_of_parent, deg_of_child;

    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;
    double running_count = 0;
    std::set<std::pair <VertexIdx,VertexIdx>> edgeList; // This list will be used to track distinct edges in the random walk.

    parent = seed;
    deg_of_parent = cg->degree(parent); // degree of parent vertex
    if (deg_of_parent == 0)  // If we are stuck with an isolate vertex, we change the seed
        parent = AlternateSeed(cg,mt);


    for (EdgeIdx wL = 0; wL < walk_length; wL++) { // Perform a random walk of length walk_length from the seed vertex
        deg_of_parent = cg->degree(parent); // degree of parent vertex
        std::uniform_int_distribution<VertexIdx> dist_parent_nbor(0, deg_of_parent - 1); // Take a step of the random walk
        EdgeIdx random_nbor_edge = cg->offsets[parent] + dist_parent_nbor(mt);
        no_of_query++;
        child = cg->nbors[random_nbor_edge]; // child is the next vertex on the random walk
        deg_of_child = cg->degree(child);// degree of the child vertex

        EdgeInfo e_struct;
        e_struct.index = random_nbor_edge;
        e_struct.src = parent;
        e_struct.dest = child;
        running_count += TriangleByEdge(cg,e_struct);

        /**
         * Update the book-keeping data structure with the visited edge and vertices
         */
        edge_list.push_back(OrderedEdge{parent, child, random_nbor_edge, deg_of_parent}); // Note that we do not care of the ordering of the edge.
        visited_edge_set[random_nbor_edge] = true;
        visited_vertex_set[parent] = true;

        /**
         * Now for the vertex with lesser degree, we visit all of its neighbors
         * and all the incident edges on it we copy the lower degree vertex in src.
         * Update relevant data structure and also update query complexity
         */
        VertexIdx src,dst;
        if (cg->degree(parent) > cg-> degree(child)) {
            src = child;
            dst = parent;
        }
        else {
            src = parent;
            dst = child;
        }
        for (EdgeIdx idx = cg->offsets[src]; idx <cg->offsets[src+1];idx++ ) {
            visited_vertex_set[cg->nbors[idx]] = true;
            visited_edge_set[idx] = true;
            EdgeIdx e2 = cg->getEdgeBinary(cg->nbors[idx],dst);
            if (e2!= -1) {
                visited_edge_set[e2] = true;
            }
            no_of_query = no_of_query + 2;
        }
        parent = child;         // The random walk proceeds with the vertex child
    }

    double edge_estimate = 0.0;
    /**
     * Depending on whether the total number of edges in the graphs is available or not, compute it.
     */
    if (!params.normalization_count_available) {
        int skip = 25;
        OrderedEdgeCollection randomEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
        Estimates edge_count_output = EstimateEdgeCount(cg, randomEdgeCollection, params, skip);
        edge_estimate = edge_count_output.estimate;
    }
    else
        edge_estimate = cg->nEdges *1.0;


    double scale = 1.0* edge_estimate / (6.0*walk_length);
    double triangleEstimate = running_count * scale ;
    //printf("p=%lf,  scale = %lf \n",p,scale);
//    printf ("Multi: %lld,  %lf\n",visited_edge_set.size(),triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;

    VertexIdx edges_seen = std::count(visited_edge_set.begin(),visited_edge_set.end(),true);
    VertexIdx vertices_seen = std::count(visited_vertex_set.begin(),visited_vertex_set.end(), true);

    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / n;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / m;
    return_estimate.query_complexity = no_of_query * 100.0 /m;

    //    return_estimate.fraction_of_vertices_seen = visited_with_nbor_vertex_set.size() * 100.0 / n;
//    return_estimate.fraction_of_edges_seen = visited_with_nbor_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}

#endif //SUBGRAPHCOUNT_SERWC_H
