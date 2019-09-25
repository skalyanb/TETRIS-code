#include <cstdlib>
#include <chrono>
#include <vector>

#include "include/GraphIO.h"
#include "include/Graph.h"
#include "include/Digraph.h"
#include "include/GetAllCounts.h"
#include "include/TriangleEstimators.h"
#include  "include/Triadic.h"
#include "include/EstimatorUtil.h"
#include "include/BaselineEstimators.h"

using namespace Escape;

// Usage:
// ./SubgraphCount input_file_path out_directory no_of_repeats seed_count sparsification_prob algo_code
// input_file_path: relative path to the input file in escape format
// out_direcory: path to the output directory . The output filename is created by the algorithm
// no_of_repeats: no of times the algo will be repeated
// seed_count: the number of start vertices for the rndom walk based algorithms
// sparsification_prob: how much of the graph we want to see? Length of the random walk
//                      is  set as sum_of_degree*sparification_prob.
// algo_code: 1 for our algo, 2 for EstTriByRWAndCount....

int main(int argc, char *argv[]) {

    Graph g;
    // Default Parameter Settings
//    int noOfRepeat = 10;
//    VertexIdx seedCount = 1;
//    double sparsification_prob = 0.1;
//    int algo_code = 1;
//
//    if (argc==1) // No config file provided, execute with default paramaters
//    {
//        std::cout << "Usage Options:\n\n";
//        std::cout << "1. ./SubgraphCount input_file_path out_directory\n\n";
//        std::cout << "2. ./SubgraphCount input_file_path out_directory no_of_repeats seed_count "
//                     "sparsification_prob algo_code \n\n";
//        std::cout << "(Option 1 runs with default parameters.)\n\n";
//    }
//    if (argc ==3 ) // Usage: ./SubgraphCount input_file_path out_directory
//    {
//        if (loadGraph(argv[1], g, 1, IOFormat::escape))
//            exit(1);
//    }
//    else if (argc == 7) // command line parameters are passed
//    {
//        if (loadGraph(argv[1], g, 1, IOFormat::escape))
//            exit(1);
//
//        printf("#Vertices = %lld, #Edges = %lld\n", g.nVertices, g.nEdges);
//
//        printf("Loaded graph from %s\n", argv[1]);
//        CGraph cg = makeCSR(g);
//        cg.sortById();
//        printf("Converted to CSR\n");
//
//        std::string filename(argv[1]);
//
//        noOfRepeat = std::strtol(argv[2],NULL,10);
//
//        EdgeIdx walkLength = g.nEdges * sparsification_prob; // g.nEdges is twice the number of edges
//        EdgeIdx subSampleSize = walkLength / 20;
//        Parameters params = {filename, seedCount, walkLength, subSampleSize, noOfRepeat, sparsification_prob};
//    }
//

    std::vector<std::string> graph_path, graph_out_dir;
    std::vector<Count> trueTriangleCount;

//    graph_path.push_back ("graphs/graph_in_edges_format/small-test.edges");
//    trueTriangleCount.push_back(8); // flickr

//    graph_path.push_back ("graphs/graph_in_edges_format/soc-flickr.edges");
//    trueTriangleCount.push_back(58771288); // flickr

    graph_path.push_back ("graphs/graph_in_edges_format/socfb-A-anon.edges");
    trueTriangleCount.push_back(55606428); //socfb-A- anon

    graph_path.push_back ("graphs/graph_in_edges_format/soc-livejournal.edges");
    trueTriangleCount.push_back(83552703); //livejournal

    graph_path.push_back ("graphs/graph_in_edges_format/soc-flickr-und.edges");
    trueTriangleCount.push_back(548658705); //flick-und

    graph_path.push_back("graphs/graph_in_edges_format/soc-orkut.edges");
    trueTriangleCount.push_back(524643952);  // orkut
//
    // The true triangle count: get by running exact count algorithm
    //soc-orkut ==524,643,952; soc-flicks == 58,771,28; soc-flickr-und = 548,658,705; livejournal = 83,552,703;
    // socfb-A-anon =  55,606,428

//    if (loadGraph(argv[1], g, 1, IOFormat::escape))
//        exit(1);
    for (int i = 0; i < graph_path.size(); i++) {

        if (loadGraph(graph_path[i].c_str(), g, 1, IOFormat::escape))
            exit(1);

        printf("#Vertices = %lld, #Edges = %lld\n", g.nVertices, g.nEdges);

        printf("Loaded graph from %s\n", graph_path[i].c_str());
        CGraph cg = makeCSR(g);
        cg.sortById();
        printf("Converted to CSR\n");

        std::string filename(graph_path[i]);
        int noOfRepeat = 5;

        std::vector <double> sparsification_prob_list {0.001,0.005,0.01,0.05,0.1};
        //std::vector <double> sparsification_prob_list {0.1,0.2,0.3,0.4};

        for ( auto sparsification_prob : sparsification_prob_list) {
            // Set up input parameters for weighted sampling estimators
            VertexIdx seedCount = 1;
            EdgeIdx walkLength = g.nEdges * sparsification_prob; // g.nEdges is twice the number of edges
            EdgeIdx subSampleSize = walkLength / 20;
            Parameters params = {filename, seedCount, walkLength, subSampleSize, noOfRepeat, sparsification_prob};

            //TriangleEstimator(&cg, params_1, params_2, trueTriangleCount[i]);

//            params.algo_name = "EstTriByRWWfgtdSamp";
//            TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriByRWandWghtedSampling);

            //TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriBySparsification);

            //TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriByUniformSampling);

            //params.algo_name = "EstTriByEdgeSampleAndCount";
            //TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriByEdgeSampleAndCount);

            //params.algo_name = "EstTriByRW";
            //TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriByRW);

            params.algo_name = "EstTriByRWAndCount";
            TriangleEstimator(&cg, params, trueTriangleCount[i], EstTriByRWAndCount);

            //CountExactTriangles (&cg);
        }
    }
    return 0;
}