//
// Created by Suman Kalyan Bera on 2020-01-26.
//

#include <cstdlib>
#include <chrono>
#include <vector>
#include <math.h>

#include "include/GraphIO.h"
#include "include/Graph.h"
#include "include/Digraph.h"
#include "include/GetAllCounts.h"
#include "include/TriangleEstimators.h"
#include  "include/Triadic.h"
#include "include/EstimatorUtil.h"
//#include "include/BaselineEstimators.h"
#include "include/util/ConfigReader.h"
#include "include/EstimateEdgeCount.h"

#include "include/baseline/VertexMCMC.h"
#include "include/baseline/SubgraohRandomWalk_SRW.h"
#include "include/baseline/SERWC.h"
#include "include/TETRIS.h"

using namespace Escape;

// Usage:
// ./SubgraphCount input_file_path no_of_repeats sparsification_prob algo_code
// input_file_path: relative path to the input file in escape format
// out_direcory: path to the output directory . The output filename is created by the algorithm
// no_of_repeats: no of times the algo will be repeated
// seed_count: the number of start vertices for the rndom walk based algorithms
// sparsification_prob: how much of the graph we want to see? Length of the random walk
//                      is  set as sum_of_degree*sparification_prob.
// algo_code: 1 for our algo, 2 for SERWC....

// The true triangle count: get by running exact count algorithm
//soc-orkut ==524,643,952; soc-flicks == 58,771,288; soc-flickr-und = 548,658,705; livejournal = 83,552,703;
// socfb-A-anon =  55,606,428 soc-sinaweibo.edges = 212,977,684
// soc-friendster.edges = 4,173,724,142
// soc-twitter-konect.edges =  34,824,916,864

int main(int argc, char *argv[]) {

    CGraph cg;
    if (argc != 2) {
        std::cout << "Usage File: ./SubgraphCount script_file\n\n";
        return 0;
    }
    config_params cfp = LoadConfig(argv[1]);

    for (int i = 0; i < cfp.input_files.size(); i++) {
        /**
         * Load the input graph and print relevant parameters
         */
        auto startTime = std::chrono::high_resolution_clock::now();
        if (loadGraphCSR(cfp.input_files[i].c_str(), cg, 1))
            exit(1);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "Time to load graph = " <<
                  std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()
                  << " seconds" << std::endl;
        std::cout << "#Graph name : " << cfp.input_files[i] << "#Vertices = " << cg.nVertices << ",#Edges = " << cg.nEdges << std::endl;

        /**
         * Load relevant parameters
         */
        Parameters params;
        params.filename = cfp.input_files[i];
        params.no_of_repeat = cfp.no_of_repeats;
        params.print_to_console = cfp.print_to_console;
        params.print_to_file = cfp.print_to_file;
        params.out_directory = cfp.out_directory;
        params.normalization_count_available = cfp.edge_count_available;

        /**
         * Set up randomness initiator
         */
        std::random_device rd;
        std::mt19937 mt(rd());

        /**
         * For each sprisification parameter and seed count and algo_name,
         * run an instance
         */
        for (auto sparsification_prob : cfp.sparsification_prob) {
            for (auto seed_count: cfp.seed_count) {
                for (auto algo_name : cfp.algo_names) {
                    /**
                     * Update parameters for this particular run
                     */
                    params.algo_name = algo_name;
                    params.sparsification_prob = sparsification_prob;
                    params.walk_length = cg.nEdges * sparsification_prob; // g.nEdges is twice the number of edges
                    params.subsample_size = params.walk_length * cfp.subsample_prob;
                    EdgeIdx triangle_count = cfp.triangle_count[i];
                    params.seed_count = seed_count;
                    params.seed_vertices.clear();

                    /**
                     * We fix the seed vertices and run for params.no_of_repeat many iterations with the
                     * seed vertex remaining fixed.
                     * We just need to
                     * sample uniform random seed vertex from the entire graph.
                     * THIS IS THE NORMAL MODE IN WHICH ALL BASELINE AND OUR
                     * ALGORITHMS ARE EXECUTED
                     */
                    VertexIdx n = cg.nVertices;
                    std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);
                    for (VertexIdx sC = 0; sC < params.seed_count; sC++) {
                        VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
                        params.seed_vertices.emplace_back(seed);
                    }

                    if (algo_name == "TETRIS") {                     // Our Algorithm
                        TriangleEstimator(&cg, params, triangle_count, TETRIS);
                    } else if (algo_name == "VertexMCMC") {
                        // Vertex MCMC sees about 2.5 time more of the edges for same walk length. For comparsion, we adjust the walk length accordingly.
                        params.walk_length = floor(params.walk_length *22 / 80);
                        TriangleEstimator(&cg, params, triangle_count, VertexMCMC);
                    }
                    else if (algo_name == "SRW1") {
                        params.walk_length = floor(params.walk_length * 22 / 40);
                        params.CSS = cfp.CSS;
                        params.NB = cfp.NB;
                        TriangleEstimator(&cg, params, triangle_count, SRW1);
                    }
                    else if (algo_name == "SERWC") { // Baseline:  do a random walk and count the number of triangles incident on each edge. Then scale.
                        params.walk_length = floor(params.walk_length * 1 / 500);
                        TriangleEstimator(&cg, params, triangle_count, SERWC);
                    }
                    else
                        std::cout << "Unknown algorithm option. \n";
//                            // Baseline: sample an edge and count the number of triangles incident on it. Then scale.
//                        else if (algo_name == "EstTriByEdgeSampleAndCount")
//                            TriangleEstimator(&cg, params, triangle_count, EstTriByEdgeSampleAndCount);
//                            // Baseline: do a random walk and count the teiangles in induces multi-graph. Scale.
//                        else if (algo_name == "EstTriByRW")
//                            TriangleEstimator(&cg, params, triangle_count, EstTriByRW);
//                            // Sample an edge, sample a neighbor and estimate the teiangles incident on the neighbor. Scale up.
//                        else if (algo_name == "EstTriByRWandNborSampling")
//                            TriangleEstimator(&cg, params, triangle_count, EstTriByRWandNborSampling);
//                            // Sample each edge with probability p and count traingles in the subsampled graph. Scale by 1/p^3.
//                        else if (algo_name == "EstTriBySparsification")
//                            TriangleEstimator(&cg, params, triangle_count, EstTriBySparsification);
//                            // Uniformly Sample edges, and count the number of triangles in the multi-graph.
//                        else if (algo_name == "EstTriByUniformSampling")
//                            TriangleEstimator(&cg, params, triangle_count, EstTriByUniformSampling);
                }
            }

        }
    }
    return 0;
}