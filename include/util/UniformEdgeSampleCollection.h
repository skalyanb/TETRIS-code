//
// Created by Suman Kalyan Bera on 2020-01-28.
//

#ifndef SUBGRAPHCOUNT_UNIFORMEDGESAMPLECOLLECTION_H
#define SUBGRAPHCOUNT_UNIFORMEDGESAMPLECOLLECTION_H

#include <random>
#include "../EstimatorUtilStruct.h"

OrderedEdgeCollection GetEdgesByUniSampling(CGraph *cg, const Parameters &params, std::mt19937 mt) {
    /**
      * Setting of various input parameters for executing the algorithm
    */
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges;
    VertexIdx seed_count = params.seed_count;
    EdgeIdx walk_length = params.walk_length;

    /**
     * Data structures for book-keeping and local variables
    */
    std::vector<OrderedEdge> edge_list;
    std::vector<bool> visited_edge_set (m,false);
    std::vector<bool> visited_vertex_set (n,false);
    VertexIdx src = 0, dst = 0;
    EdgeIdx nEdges = 0;
    /**
       * Set up random number generator
    */
    std::uniform_int_distribution<EdgeIdx> unif_rand_edge(0, m - 1);

    /**
     * Sample uniform random edges and store them in the structure
     */
    for (EdgeIdx i = 0; i < walk_length ; i++) {
        EdgeIdx e = unif_rand_edge(mt);
        dst = cg->nbors[e];
        // binary search the array cg->offset to find the index v, such that edges[e] lies between cg->offset[v]
        // and cg->offset[v+1]. We achieve this by calling std::upper_bound
        src = std::upper_bound (cg->offsets, cg->offsets+n, e) - cg->offsets -1;

        // sanity check: the edge {src,dst} must have index edges[e].
        if (cg->getEdgeBinary(src,dst) != e)
            printf("Bug in the code!! Run!!!! \n\n");

        /**
          * Update the edge collection data structure with relevant information about the edge
         */
        VertexIdx u, v, low_deg;
        if (cg->degree(dst) < cg->degree(src) || (cg->degree(dst) == cg->degree(src) && dst < src)) {
            u = dst;
            v = src;
            low_deg = cg->degree(dst);
        } else {
            u = src;
            v = dst;
            low_deg = cg->degree(src);
        }
        edge_list.push_back(OrderedEdge{u, v, e, low_deg});

        /**
         * Update the book-keeping data structure with the visited edge and vertices
         */
        visited_edge_set[e] = true;
        visited_vertex_set[src] = true;
        visited_vertex_set[dst] = true;
    }

    OrderedEdgeCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
    return returnEdgeCollection;
}

#endif //SUBGRAPHCOUNT_UNIFORMEDGESAMPLECOLLECTION_H
