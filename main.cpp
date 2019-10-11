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
#include "include/BaselineEstimators.h"
#include "include/ConfigReader.h"
#include "include/EstimateEdgeCount.h"

using namespace Escape;

// Usage:
// ./SubgraphCount input_file_path no_of_repeats sparsification_prob algo_code
// input_file_path: relative path to the input file in escape format
// out_direcory: path to the output directory . The output filename is created by the algorithm
// no_of_repeats: no of times the algo will be repeated
// seed_count: the number of start vertices for the rndom walk based algorithms
// sparsification_prob: how much of the graph we want to see? Length of the random walk
//                      is  set as sum_of_degree*sparification_prob.
// algo_code: 1 for our algo, 2 for EstTriByRWAndCountPerEdge....

// The true triangle count: get by running exact count algorithm
//soc-orkut ==524,643,952; soc-flicks == 58,771,288; soc-flickr-und = 548,658,705; livejournal = 83,552,703;
// socfb-A-anon =  55,606,428 soc-sinaweibo.edges = 212,977,684
// soc-friendster.edges = 4,173,724,142
// soc-twitter-konect.edges =  34,824,916,864

int main(int argc, char *argv[]) {

    Graph g;
    if (argc != 2) {
        std::cout << "Usage File: ./SubgraphCount script_file\n\n";
        return 0;
    }
    config_params cfp = LoadConfig(argv[1]);

    for (int i = 0; i < cfp.input_files.size(); i++) {

        //Uplaod the graph from input path
        if (loadGraph(cfp.input_files[i].c_str(), g, 1, IOFormat::escape))
            exit(1);

        printf("#Vertices = %lld, #Edges = %lld\n", g.nVertices, g.nEdges);

        printf("Loaded graph from %s\n", cfp.input_files[i].c_str());
        CGraph cg = makeCSR(g);
        cg.sortById();
        printf("Converted to CSR\n");

        Parameters params;
        params.filename = cfp.input_files[i];
        params.no_of_repeat = cfp.no_of_repeats;
        params.print_to_console = cfp.print_to_console;
        params.print_to_file = cfp.print_to_file;

        // For each sprisification parameter and seed count and algo_name, run an instance
        for (auto algo_name : cfp.algo_names) {
            params.algo_name = algo_name;

            for (auto sparsification_prob : cfp.sparsification_prob) {
                params.sparsification_prob = sparsification_prob;
                params.walk_length = g.nEdges * sparsification_prob; // g.nEdges is twice the number of edges
                params.subsample_size = params.walk_length * cfp.subsample_prob;
                EdgeIdx triangle_count = cfp.triangle_count[i];

                for (auto seed_count: cfp.seed_count) {
                    params.seed_count = seed_count;
                    params.seed_vertices.clear();
                    // We fix the seed vertices and for run params.no_of_repeat many iterations with the
                    // seed vertex remaining fixed.
                    std::random_device rd;
                    std::mt19937 mt(rd());
                    // If degree_bin_seed is false, then we just need to
                    // sample uniform random seed vertex from the entire graph.
                    // THIS IS THE NORMAL MODE IN WHICH ALL BASELINE AND OUR
                    // ALGORITHMS ARE EXECUTED
                    if (cfp.degree_bin_seed == false) {
                        VertexIdx n = cg.nVertices;
                        std::uniform_int_distribution<VertexIdx> dist_seed_vertex(0, n - 1);
                        for (VertexIdx sC = 0; sC < params.seed_count; sC++) {
                            VertexIdx seed = dist_seed_vertex(mt); // TODO: verify randomness
                            params.seed_vertices.emplace_back(seed);
                        }
                        // Our Algorithm
                        if (algo_name == "EstTriByRWandWghtedSampling"){
                            params.algo_name = "_new_" + algo_name;
                            TriangleEstimator(&cg, params, triangle_count, EstTriByRWandWghtedSampling);
                        }
                            // Baseline: sample an edge and count the number of triangles incident on it. Then scale.
                        else if (algo_name == "EstTriByEdgeSampleAndCount")
                            TriangleEstimator(&cg, params, triangle_count, EstTriByEdgeSampleAndCount);
                            // Baseline:  do a random walk and count the number of triangles incident on each edge. Then scale.
                        else if (algo_name == "EstTriByRWAndCountPerEdge")
                            TriangleEstimator(&cg, params, triangle_count, EstTriByRWAndCountPerEdge);
                            // Baseline: do a random walk and count the teiangles in induces multi-graph. Scale.
                        else if (algo_name == "EstTriByRW")
                            TriangleEstimator(&cg, params, triangle_count, EstTriByRW);
                            // Sample an edge, sample a neighbor and estimate the teiangles incident on the neighbor. Scale up.
                        else if (algo_name == "EstTriByRWandNborSampling")
                            TriangleEstimator(&cg, params, triangle_count, EstTriByRWandNborSampling);
                            // Sample each edge with probability p and count traingles in the subsampled graph. Scale by 1/p^3.
                        else if (algo_name == "EstTriBySparsification")
                            TriangleEstimator(&cg, params, triangle_count, EstTriBySparsification);
                            // Uniformly Sample edges, and count the number of triangles in the multi-graph.
                        else if (algo_name == "EstTriByUniformSampling")
                            TriangleEstimator(&cg, params, triangle_count, EstTriByUniformSampling);
                        else if (algo_name == "EdgeEstimator")
                            EdgeEstimator(&cg, params);
                        else
                            std::cout << "Unknown algorithm option. \n";
                    }
                    // Otherwise, we iterate over all the vertices, and sample 5 vertex from each degree range
                    // to figure out the degree range, we take the log of the degree.
                    // WE ONLY NEED TO DO THIS FOR OUR ALGORITHM.
                    else {
                        params.seed_vertices.emplace_back(-1); // dummy place holder for now
                        VertexIdx num_of_bucket = floor(log10(cg.nVertices));
                        std::vector<std::vector<VertexIdx >> deg_bins(num_of_bucket);
                        std::vector<std::vector<VertexIdx >> random_seeds(num_of_bucket);
                        for (VertexIdx v = 0; v < cg.nVertices; v++) {
                            VertexIdx deg = cg.degree(v);
                            if (deg != 0) {
                                VertexIdx bucket_id = floor(log10(deg));
                                deg_bins[bucket_id].emplace_back(v);
                            }
                        }
                        // Now sample 5 vertices uniformly at random from each bin (with replacement)
                        for (int nb = 0; nb < num_of_bucket; nb++) {
                            if (!deg_bins[nb].empty()) {
                                std::uniform_int_distribution<VertexIdx> dist_bucket_seed_vertex(0, deg_bins[nb].size() - 1);
                                for (int j = 0; j < 5; j++) {
                                    VertexIdx random_seed = dist_bucket_seed_vertex(mt);
                                    random_seeds[nb].emplace_back(deg_bins[nb][random_seed]);
                                }
                            }
                        }
                        // Clear the excess memory created by this seeding process.
                        std::vector<std::vector<VertexIdx >>().swap(deg_bins);
                        if (algo_name == "EstTriByRWandWghtedSampling"){
                            for (int nb = 0; nb < num_of_bucket; nb++) {
                                if (!random_seeds[nb].empty()) {
                                    for (int j = 0; j < 4; j++) {
                                        VertexIdx seed = random_seeds[nb][j];
                                        params.seed_vertices[0]= seed;
                                        params.algo_name = "EstTriByRWandWghtedSampling_" + std::to_string(nb);
                                        TriangleEstimator(&cg, params, triangle_count, EstTriByRWandWghtedSampling);
                                    }
                                }
                            }
                        }
                        else
                            std::cout << "Unknown algorithm option. \n";
                    }
                }
            }

        }
    }
    return 0;
}