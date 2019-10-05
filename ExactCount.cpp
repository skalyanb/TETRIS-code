//
// Created by Suman Kalyan Bera on 2019-09-28.
//

#include <iostream>

#include "include/GraphIO.h"
#include "include/Graph.h"
#include "include/TriangleEstimators.h"

using namespace Escape;
int main(int argc, char *argv[]) {

    if (argc!=2) {
        std::cout << "Usage: ./ExactCount input_file_name";
        return -1;
    }

    //Uplaod the graph from input path
    Graph g;
    if (loadGraph(argv[1], g, 1, IOFormat::escape))
        exit(1);
    printf("#Vertices = %lld, #Edges = %lld\n", g.nVertices, g.nEdges);

    printf("Loaded graph from %s\n", argv[1]);
    CGraph cg = makeCSR(g);
    cg.sortById();
    printf("Converted to CSR\n");

    Estimates out = CountExactTriangles(&cg);
    std::cout << "#Vertice =" << cg.nVertices << " #Edges =" << cg.nEdges << " #Triangles ="  << out.triangle_estimate  <<"\n";
    return 0;
}