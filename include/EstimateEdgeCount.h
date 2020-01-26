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


OrderedEdgeCollection GetEdgesByUniSampling(CGraph *cg, Parameters params, std::mt19937 mt) {
    // Get the number of vertices, number of edges, Initialize scaling to number of edges.
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    VertexIdx seed_count = params.seed_count;
    EdgeIdx walk_length = params.walk_length;

    // The list of edges found by the random walk process
    std::vector<OrderedEdge> edge_list;

    // Keep track of the vertices and edges seen so far by the random walk
    std::vector<bool> visited_edge_set (m,false);
    std::vector<bool> visited_vertex_set (n,false);

    VertexIdx src = 0, dst = 0;
    EdgeIdx nEdges = 0;
    std::uniform_int_distribution<EdgeIdx> unif_rand_edge(0, m - 1);


    // Initialize a uniform distribution for generating seed vertices
    //std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);

    // Perform a random walk for seed_count many times
    for (EdgeIdx i = 0; i < walk_length ; i++) {
        EdgeIdx e = unif_rand_edge(mt);
        dst = cg->nbors[e];
        // binary search the array cg->offset to find the index v, such that edges[e] lies between cg->offset[v]
        // and cg->offset[v+1]. We achieve this by calling std::upper_bound
        src = std::upper_bound (cg->offsets, cg->offsets+n, e) - cg->offsets -1;

        // sanity check: the edge {src,dst} must have index edges[e].
        if (cg->getEdgeBinary(src,dst) != e)
            printf("Bug in the code!! Run!!!! \n\n");

        edge_list.push_back(OrderedEdge{src, dst, e, 0});
        visited_edge_set[e] = true;
        visited_vertex_set[src] = true;
        visited_vertex_set[dst] = true;
    }

    // Create return structure for this function.
    // The number of edges produced by the random walk is walk_length;
//    std::copy(visited_edge_list.begin(),
//              visited_edge_list.end(),
//              std::inserter(visited_edge_set, visited_edge_set.end()));
//
//    std::copy(visited_vertex_list.begin(),
//              visited_vertex_list.end(),
//              std::inserter(visited_vertex_set, visited_vertex_set.end()));

    OrderedEdgeCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
//    printf("Random walk: walk length = %lld ",walk_length);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",visited_edge_set.size(),visited_vertex_set.size());

    return returnEdgeCollection;
}


Estimates EstimateEdgeCount (CGraph *cg, OrderedEdgeCollection randomEdgeCollection, Parameters params, int skip)
{
//    OrderedEdgeCollection randomEdgeCollection;

    // Set up random number generator
//    std::random_device rd;
//    // Using this random number generator initializize a PRNG: this PRNG is passed along to
//    // draw an element from various distribution
//    std::mt19937 mt(rd());

    EdgeIdx collsion_count =0;
    double numerator = 0;
    int c = 1;

    // Now count the number of collision in the random edge collection
    //1. By random walk
//    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    c = skip;
    //2. By uniform edge collection
//    randomEdgeCollection = GetEdgesByUniSampling(cg, params, mt);

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

    printf("True edge count = %lld, initial sample size = %lld\n",cg->nEdges,randomEdgeCollection.edge_list.size());
    printf("Num edges=%lld, numerator=%lf, numerator2=%lf, collision=%lld,edge_estimate=%lf.\n",
            num_edge, numerator,numerator, collsion_count,global_estimate);
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
