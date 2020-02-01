//
// Created by Suman Kalyan Bera on 2020-01-31.
//

#ifndef SUBGRAPHCOUNT_UESS_H
#define SUBGRAPHCOUNT_UESS_H

#include <random>
#include "../EstimatorUtilStruct.h"
#include "../util/BaselineUtil.h"
#include "../util/UniformEdgeSampleCollection.h"

/**
 * Algorithm:
 * 1. Sample R many edges from the graph, independently and uniformly at random.
 * 2. Let the subsampled graph be G_R.
 * 3. Count the number of triangles in G_R and scale the count up by 1/p^3 where p = R/m.
 * @param cg
 * @param params
 * @return
 *
 * This algorithm is implemented by a three step process.
 * First, we construct a new graph object from the subsampled edges; call it G_p
 * Then, we convert it into a CGraph object
 * Finally, we call exact triangle count routine on CGraph object
 */

Estimates UESS(CGraph *cg, Parameters params)
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

    std::vector<VertexIdx > srcs;
    std::vector<VertexIdx > dsts;

    for (auto edge : randomEdgeCollection.edge_list) {
        srcs.emplace_back(edge.u);
        dsts.emplace_back(edge.v);
        srcs.emplace_back(edge.v);  // The graph representation (adjacency list) requires both the edges to be present
        dsts.emplace_back(edge.u);
    }
    /**
     * Construct the induced graph G_p
     * Note that G_p is a multi-graph. When counting triangles exactly in G_p, we consider this multiplicity.
     */
    CGraph CG_p = MakeMultiGraph(cg->nVertices,srcs,dsts);

    /**
     * Count the exact number of triangles in the sparsified graph CG_p and scale it up.
     */

    Estimates triangles_in_Gp = CountExactTriangles(&CG_p);
    delCGraph(CG_p);
    double p = randomEdgeCollection.no_of_edges *2.0 / cg->nEdges;
    double scale = 1.0 / pow(p,3);
    double triangleEstimate = triangles_in_Gp.estimate * scale ;

    Estimates return_estimate = {};
    return_estimate.estimate = triangleEstimate;

    VertexIdx edges_seen = std::count(randomEdgeCollection.visited_edge_set.begin(),randomEdgeCollection.visited_edge_set.end(),true);
    VertexIdx vertices_seen = std::count(randomEdgeCollection.visited_vertex_set.begin(),randomEdgeCollection.visited_vertex_set.end(), true);

    return_estimate.fraction_of_vertices_seen = vertices_seen * 100.0 / cg->nVertices;
    return_estimate.fraction_of_edges_seen = edges_seen * 100.0 / cg->nEdges;
    return_estimate.query_complexity = randomEdgeCollection.no_of_query * 100.0 /cg->nEdges;

    return return_estimate;
}


#endif //SUBGRAPHCOUNT_UESS_H
