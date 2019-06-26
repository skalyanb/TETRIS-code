#include <cstdlib>
#include <chrono>
#include <vector>

#include "include/GraphIO.h"
#include "include/Graph.h"
#include "include/Digraph.h"
#include "include/GetAllCounts.h"
#include "include/Triangle.h"

using namespace Escape;

/**
 *
 * @param argc
 * @param argv
 * @return
 */

int main(int argc, char *argv[])
{
    Graph g;
    if (loadGraph(argv[1], g, 1, IOFormat::escape))
        exit(1);

    printf("#Vertices = %lld, #Edges = %lld\n", g.nVertices, g.nEdges/2);

    printf("Loaded graph\n");
    CGraph cg = makeCSR(g);
    cg.sortById();
    printf("Converted to CSR\n");

    // Set up input parameters for estimators
    VertexIdx seedCount = 1;
    EdgeIdx walkLength = g.nEdges/100; // g.nEdges is twice the number of edges
    EdgeIdx subSampleSize = walkLength/50;
    int noOfRepeat = 100;
    std::string filename(argv[1]);
    Parameters params = {seedCount,walkLength,subSampleSize,noOfRepeat,filename};

    // The true triangle count: get by running exact count algorithm
    //soc-orkut ==524,643,952
    //soc-flicks == 58,771,288
    // soc-flickr-und = 548,658,705
    //Count trueTriangleCount = 548658705;
    Count trueTriangleCount = 524643952;
    //Count trueTriangleCount = 58771288;


    TriangleEstimator (&cg, params, trueTriangleCount);

    return 0;
}