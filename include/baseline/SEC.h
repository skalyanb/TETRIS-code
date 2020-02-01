//
// Created by Suman Kalyan Bera on 2020-01-31.
//

#ifndef SUBGRAPHCOUNT_SEC_H
#define SUBGRAPHCOUNT_SEC_H


#include <set>
#include <unordered_map>

#include "../EstimatorUtilStruct.h"
#include "../TriangleEstimators.h"
#include "../util/UniformEdgeSampleCollection.h"
#include "../util/BaselineUtil.h"
#include "../EstimateEdgeCount.h"

#include "../util/RandomWalkUtils.h"

/**
 * Algorithm Details:
 * 1. Sample an edge uniformly at random from the graph
 * 2. Count the number of triangles incident on the edge.
 * 3. Repeat the process and take average count.
 * 4. Scale up by m/3 where m is the number of edges.
 * @param cg
 * @param params
 * @return
 */

Estimates SEC(CGraph *cg, Parameters params)
{
    /**
    * Set up local variables. skip parameter is used in estimating the number of edges in the graph.
    */
    OrderedEdgeCollection randomEdgeCollection;
    double running_count = 0.0, edge_estimate=0.0;
    int skip = 1;

    /**
    * Set up random number generator
    * Using this random number generator initializize a PRNG: this PRNG is passed along to
    * draw an element from various distributions
    */
    std::random_device rd;
    std::mt19937 mt(rd());

    /**
     * Perform a random walk through the graph, and collect the edges in the process
     * The fourth parameter is set to false to indicate that we do not care about the
     * vertex ordering of each edge
     */
    randomEdgeCollection = GetEdgesByUniSampling(cg, params, mt);

    /**
     * Depending on whether the total number of edges in the graphs is available or not, compute it.
     */
    if (!params.normalization_count_available) {
        Estimates edge_count_output = EstimateEdgeCount(cg, randomEdgeCollection, params, skip);
        edge_estimate = edge_count_output.estimate;
    }
    else
        edge_estimate = cg->nEdges *1.0;

    /** We now go over the edges collected by the random walk and count the
     * number of triangles incident on each edge.
     * However, we can not do this for every edge as it would be too expensive.
     * Instead, we perform this task for a small fraction of the total edges in the random walk
     */
    int sparsity = 500, no_of_samples=0;
    EdgeIdx total_edge = randomEdgeCollection.edge_list.size();

    for (VertexIdx i=0; i<total_edge; i=i+sparsity) {
        EdgeInfo e_struct;
        VertexIdx parent,child;
        parent = randomEdgeCollection.edge_list[i].u;
        child = randomEdgeCollection.edge_list[i].v;

        e_struct.index = randomEdgeCollection.edge_list[i].index;
        e_struct.src = parent;
        e_struct.dest = child;
        running_count += TriangleByEdge(cg,e_struct);
        no_of_samples++;
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
            randomEdgeCollection.visited_vertex_set[cg->nbors[idx]] = true;
            randomEdgeCollection.visited_edge_set[idx] = true;
            EdgeIdx e2 = cg->getEdgeBinary(cg->nbors[idx],dst);
            if (e2!= -1) {
                randomEdgeCollection.visited_edge_set[e2] = true;
            }
            randomEdgeCollection.no_of_query = randomEdgeCollection.no_of_query + 2;
        }
    }



    double scale = 1.0* edge_estimate / (6.0*no_of_samples);
    double triangleEstimate = running_count * scale ;
    //printf("p=%lf,  scale = %lf \n",p,scale);
//    printf ("Multi: %lld,  %lf\n",visited_edge_set.size(),triangleEstimate);

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;

    VertexIdx edges_seen = std::count(randomEdgeCollection.visited_edge_set.begin(),randomEdgeCollection.visited_edge_set.end(),true);
    VertexIdx vertices_seen = std::count(randomEdgeCollection.visited_vertex_set.begin(),randomEdgeCollection.visited_vertex_set.end(), true);

    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / cg->nEdges;
    return_estimate.query_complexity = randomEdgeCollection.no_of_query * 100.0 /cg->nEdges;

    //    return_estimate.fraction_of_vertices_seen = visited_with_nbor_vertex_set.size() * 100.0 / n;
//    return_estimate.fraction_of_edges_seen = visited_with_nbor_edge_set.size() * 100.0 / m;
    // Why the multiplication by 2? The fraction of edges seen is effectively
    // compared against all the entries in the adjacency list, which is 2m.

    return return_estimate;
}


#endif //SUBGRAPHCOUNT_SEC_H
