//
// Created by Suman Kalyan Bera on 2019-10-08.
//

#ifndef SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H
#define SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H

#include <cmath>
#include <map>
#include "TriangleEstimators.h"
#include "EstimatorUtil.h"

struct EdgeEstimates {
    double edge_estimate = 0.0;
    double fraction_of_vertices_seen = 1.0;
    double fraction_of_edges_seen = 1.0; // This fraction is with respect to twice the number of edges
};

OrderedEdgeCollection GetEdgesByUniSampling(CGraph *cg, Parameters params, std::mt19937 mt) {
    // Get the number of vertices, number of edges, Initialize scaling to number of edges.
    VertexIdx n = cg->nVertices;
    EdgeIdx m = cg->nEdges; // TODO note that m is double of the number of edges in the graph
    VertexIdx seed_count = params.seed_count;
    EdgeIdx walk_length = params.walk_length;

    // The list of edges found by the random walk process
    std::vector<OrderedEdge> edge_list;

    // Keep track of the vertices and edges seen so far by the random walk
    std::unordered_set<EdgeIdx> visited_edge_set;
    std::unordered_set<VertexIdx> visited_vertex_set;

    std::vector<EdgeIdx> visited_edge_list;
    std::vector<VertexIdx> visited_vertex_list;

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
        visited_edge_list.emplace_back(e);
        visited_vertex_list.emplace_back(src);
        visited_vertex_list.emplace_back(dst);
    }

    // Create return structure for this function.
    // The number of edges produced by the random walk is walk_length;
    std::copy(visited_edge_list.begin(),
              visited_edge_list.end(),
              std::inserter(visited_edge_set, visited_edge_set.end()));

    std::copy(visited_vertex_list.begin(),
              visited_vertex_list.end(),
              std::inserter(visited_vertex_set, visited_vertex_set.end()));

    OrderedEdgeCollection returnEdgeCollection = {walk_length, edge_list, visited_edge_set, visited_vertex_set};
//    printf("Random walk: walk length = %lld ",walk_length);
//    printf("Edges Seen=%lu, Vertices Seen=%lu\n",visited_edge_set.size(),visited_vertex_set.size());

    return returnEdgeCollection;
}



Estimates EstimateEdgeCount (CGraph *cg, Parameters params)
{
    OrderedEdgeCollection randomEdgeCollection;

    // Set up random number generator
    std::random_device rd;
    // Using this random number generator initializize a PRNG: this PRNG is passed along to
    // draw an element from various distribution
    std::mt19937 mt(rd());

    EdgeIdx collsion_count =0;
    double numerator = 0;
    int c = 1;

    // Now count the number of collision in the random edge collection
    //1. By random walk
    randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
    c = 25;
    //2. By uniform edge collection
//    randomEdgeCollection = GetEdgesByUniSampling(cg, params, mt);

    EdgeIdx num_edge = randomEdgeCollection.edge_list.size();
    std::vector<EdgeIdx > edge_list;
    // Copy every c-th element from the vector and count collision in them
    for (EdgeIdx i =0; i < num_edge; i=i+c ) {
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
        collsion_count += p.second * (p.second-1) /2;
    }



//    for (EdgeIdx i =0; i < num_edge; i++ ) {
//        EdgeIdx e_i = randomEdgeCollection.edge_list[i].index;
//        for (EdgeIdx j = i+1; j< num_edge; j++) {
//            EdgeIdx e_j = randomEdgeCollection.edge_list[j].index;
//            numerator++;
//            if ( e_i == e_j)
//                collsion_count++;
//        }
//    }
    numerator = num_edge * (num_edge-1) /2.0;
    double edge_estimate = 1.0 * numerator / collsion_count;
    printf("True edge count = %lld, initial sample size = %lld\n",cg->nEdges,randomEdgeCollection.edge_list.size());
    printf("Num edges=%lld, numerator=%lf, numerator2=%lf, collision=%lld,edge_estimate=%lf.\n",
            num_edge, numerator,numerator, collsion_count,edge_estimate);
    Estimates output;
    output.triangle_estimate = edge_estimate; // This is so bad! the triangle count name is actually stroing the edge count.
                                                // We definitely need to fix this
    output.fraction_of_vertices_seen = randomEdgeCollection.visited_vertex_set.size() * 100.0 / cg->nVertices;
    output.fraction_of_edges_seen = randomEdgeCollection.visited_edge_set.size() * 100.0 / cg->nEdges;

    return output;
}

void EdgeEstimator (CGraph *cg, Parameters params) {
    std::vector<Estimates> estimates;

    for (Count i = 0; i < params.no_of_repeat; i++) {
        Estimates est = EstimateEdgeCount(cg, params);
        estimates.push_back(est);
        std::cout << i << "\n";
    }
    EstimatorStats est_stats = GetErrorStatistics(estimates, cg->nEdges);

    // Print to console
    if (params.print_to_console) {
        WriteHeaderInOutput(stdout, params, cg, cg->nEdges);
        WriteAlgorithmOutput(stdout, params.algo_name, params, est_stats);
        WriteRawData(stdout,estimates);
    }

    // print to file
    if (params.print_to_file) {
        // take current timestamp
        std::string current_time = GetTimestamp();
        // extract the input filename
        std::string out_filename = params.filename.substr(params.filename.find_last_of("/\\") + 1);

        // create the directory output/filename/algoname, if it already does not exist
        std::string directory = "output/" + out_filename + "/" + params.algo_name;
        DIR* dir = opendir(directory.c_str());  // try to open the directory
        if (dir) {
            // Directory exists, go on to create a file in this location
            closedir(dir);
        }
        else if (ENOENT == errno){
            // Directory does not exist, try creating one
            std::string mkdir = "mkdir -p " + directory;
            if (std::system(mkdir.c_str()) == -1) {
                printf("Could not create directory output/%s/%s\n", out_filename.c_str(), params.algo_name.c_str());
                return;
            }

        }
        else {
            printf("Could not create/find directory output/%s/%s\n",out_filename.c_str(),params.algo_name.c_str());
            return;
        }

        out_filename = "output/" + out_filename + "/" + params.algo_name + "/" +
                       std::to_string(params.sparsification_prob) + "-" +
                       current_time + "-" + out_filename + "-" + params.algo_name + + ".txt";
        FILE *f = fopen(out_filename.c_str(), "w");
        if (!f) {
            printf("Could not write output. Please check for write permissions. Lcaoltion: %s\n",out_filename.c_str());
            return;
        }
        WriteHeaderInOutput(f, params, cg, cg->nEdges);
        WriteAlgorithmOutput(f, params.algo_name, params, est_stats);
        WriteRawData(f, estimates);

        fclose(f);
    }

}

#endif //SUBGRAPHCOUNT_ESTIMATEEDGECOUNT_H
