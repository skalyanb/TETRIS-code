//
// Created by Suman Kalyan Bera on 2019-10-08.
//

#ifndef SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H
#define SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H

#include <cmath>
#include <map>
//#include "TriangleEstimators.h"
//#include "EstimatorUtil.h"
#include "EstimatorUtilStruct.h"

/**
 * Algorithm Details:
 * 1. Perform a long random walk. Let the set be R. This is given as input to this method.
 * 2. Let R_skip be the set of edges taken at an interval of "skip" from the set R
 * 3.
 * @param cg
 * @param randomEdgeCollection
 * @param params
 * @param skip
 * @return
 */


Estimates EstimateEdgeCount (CGraph *cg, OrderedEdgeCollection randomEdgeCollection, Parameters params, int skip)
{
    EdgeIdx collsion_count =0;
    double numerator = 0;
    int c=skip;

    EdgeIdx total_edge = randomEdgeCollection.edge_list.size();
    EdgeIdx num_edge = 0;
    int starting_pos = 0;
    double global_estimate=0;
    for (starting_pos = 0; starting_pos < c; starting_pos++) {
        EdgeIdx local_collsion_count = 1;
        std::vector<EdgeIdx > edge_list;
        // Copy every c-th element from the vector and count collision in them
        for (EdgeIdx i =starting_pos; i < total_edge; i=i+c ) {
            EdgeIdx e_i = randomEdgeCollection.edge_list[i].index;
            edge_list.emplace_back(e_i);
        }
        num_edge = edge_list.size();
        // Find the collisons now
        //std::sort(edge_list.begin(),edge_list.end());
        std::map<EdgeIdx , VertexIdx > freq_map;
        for (auto const & x : edge_list)
            ++freq_map[x];
        for (auto const & p : freq_map) {
            // How many collisions? If frequency of an edge is k, then k(k-1)/2 many collisions
            local_collsion_count += p.second * (p.second-1) /2;
            collsion_count += local_collsion_count;
        }
        numerator = num_edge * (num_edge-1) /2.0;
        if (local_collsion_count ==0)
            printf("No!");
        double edge_estimate = 1.0 * numerator / local_collsion_count;
        global_estimate += edge_estimate;
    }
    global_estimate = 1.0 * global_estimate/c;
    collsion_count /=c;

//    printf("True edge count = %lld, initial sample size = %lld\n",cg->nEdges,randomEdgeCollection.edge_list.size());
//    printf("Num edges=%lld, numerator=%lf, numerator2=%lf, collision=%lld,edge_estimate=%lf.\n",
//            num_edge, numerator,numerator, collsion_count,global_estimate);
    Estimates output;
    output.estimate = global_estimate; // This is so bad! the triangle count name is actually stroing the edge count.
                                                // We definitely need to fix this
    VertexIdx vertices_seen = std::count(randomEdgeCollection.visited_vertex_set.begin(),randomEdgeCollection.visited_vertex_set.end(),true);
    EdgeIdx edges_seen = std::count(randomEdgeCollection.visited_edge_set.begin(),randomEdgeCollection.visited_edge_set.end(),true);
    output.fraction_of_vertices_seen = vertices_seen * 100.0 / cg->nVertices;
    output.fraction_of_edges_seen = edges_seen * 100.0 / cg->nEdges;

    return output;
}


#endif //SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H
