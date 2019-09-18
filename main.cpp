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
using namespace Escape;

int main(int argc, char *argv[]) {
    Graph g;
    std::vector<std::string> graph_path;
    std::vector<Count> trueTriangleCount;

//    CountTriangle(argc,argv);
//    return 0;
//    graph_path.push_back ("../graphs/soc-flickr.edges");
//    trueTriangleCount.push_back(58771288); // flickr
//
//    graph_path.push_back ("../graphs/soc-livejournal.edges");
//    trueTriangleCount.push_back(83552703); //livejournal
//
//    graph_path.push_back ("../graphs/soc-flickr-und.edges");
//    trueTriangleCount.push_back(548658705); //flick-und
//
//    graph_path.push_back ("../graphs/socfb-A-anon.edges");
//    trueTriangleCount.push_back(55606428); //socfb-A- anon

    graph_path.push_back("../graphs/graph_in_edges_format/soc-orkut.edges");
    trueTriangleCount.push_back(524643952);  // orkut

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
        int noOfRepeat = 20;
        // Set up input parameters for simple sampling estimators
        VertexIdx seedCount_1 = 1;
        EdgeIdx walkLength_1 = g.nEdges / 1000; // g.nEdges is twice the number of edges
        Parameters params_1 = {seedCount_1, walkLength_1, 0, noOfRepeat, filename};

        // Set up input parameters for weighted sampling estimators
        VertexIdx seedCount_2 = 1;
        EdgeIdx walkLength_2 = g.nEdges / 1000; // g.nEdges is twice the number of edges
        EdgeIdx subSampleSize = walkLength_2 / 10;
        Parameters params_2 = {seedCount_2, walkLength_2, subSampleSize, noOfRepeat, filename};

        //TriangleEstimator(&cg, params_1, params_2, trueTriangleCount[i]);
        TriangleEstimator(&cg, params_1,trueTriangleCount[i],EstTriByEdgeSparsification);
        //CountExactTriangles (&cg);
    }
    return 0;
}