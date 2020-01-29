//
// Created by Suman Kalyan Bera on 2020-01-28.
//

#ifndef SUBGRAPHCOUNT_DEGREESQUARESUM_H
#define SUBGRAPHCOUNT_DEGREESQUARESUM_H

#include <cmath>
#include <map>
#include <algorithm>

#include "EstimatorUtilStruct.h"
#include "EstimateEdgeCount.h"

/**
 * Algorithm Details:
 * This algorithm estimates \sum_v d(v) (d(v)-1)/2 from the random walk perform by VertexMCMC
 * Note the VertexMCMC utilizes a metropolis-hastings algorithm with stationary distribution
 * proportional to \binom{d(v),2}.
 * @param cg
 * @param randomEdgeCollection
 * @param params
 * @param skip
 * @return
 */


Estimates DegreeSumSquare(CGraph *cg, OrderedEdgeCollection randomEdgeCollection, Parameters params, int skip)
{
    int c = skip;

    EdgeIdx total_edge = randomEdgeCollection.edge_list.size();
    int starting_pos = 0;
    double global_estimate=0;
    std::vector <double> global_estimate_runs(c,0.0);
    for (starting_pos = 0; starting_pos < c; starting_pos++) {
        VertexIdx num_samples = 0;
        for (EdgeIdx i =starting_pos; i < total_edge; i=i+c ) {
            num_samples++;
            VertexIdx d_i = cg->degree(randomEdgeCollection.edge_list[i].u);
            global_estimate_runs[starting_pos] += (d_i -1)/2.0;
        }
        global_estimate_runs[starting_pos] = global_estimate_runs[starting_pos] * 1.0 / num_samples;
    }

    std::sort(global_estimate_runs.begin(),global_estimate_runs.end());
    global_estimate = global_estimate_runs[c/2];

    Estimates edge_count_output = EstimateEdgeCount(cg, randomEdgeCollection, params, skip);
    double edge_estimate = edge_count_output.estimate;
    global_estimate = 1.0 * edge_estimate * global_estimate;

    Estimates output;
    output.estimate = global_estimate;
    VertexIdx vertices_seen = std::count(randomEdgeCollection.visited_vertex_set.begin(),randomEdgeCollection.visited_vertex_set.end(),true);
    EdgeIdx edges_seen = std::count(randomEdgeCollection.visited_edge_set.begin(),randomEdgeCollection.visited_edge_set.end(),true);
    output.fraction_of_vertices_seen = vertices_seen * 100.0 / cg->nVertices;
    output.fraction_of_edges_seen = edges_seen * 100.0 / cg->nEdges;

//    double wedge_count = 0;
//    for (VertexIdx src=0; src < cg->nVertices; src++) {
//        VertexIdx deg = cg->degree(src);
//        wedge_count += deg * (deg-1) * 0.5;
//    }
//    double global_estimate_with_m = global_estimate * cg->nEdges / (1.0 * edge_estimate);

//    printf("True Estimate : %.3f, Estimate with Edges known : %.3fd, Final Estimate :%.3f \n", wedge_count,
//              global_estimate_with_m,
//              global_estimate);
    return output;
}
#endif //SUBGRAPHCOUNT_DEGREESQUARESUM_H
