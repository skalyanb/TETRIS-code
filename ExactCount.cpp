//
// Created by Suman Kalyan Bera on 2019-09-28.
//

#include <iostream>

#include "include/GraphIO.h"
#include "include/Graph.h"
#include "include/TriangleEstimators.h"
#include "include/EstimatorUtil.h"

using namespace Escape;

void write_output (std::string filename, CGraph *cg, Count triangle_count) {

    // take current timestamp
    std::string current_time = GetTimestamp();
    // extract the input filename
    std::string out_filename = filename.substr(filename.find_last_of("/\\") + 1);

    // create the directory output/filename/algoname, if it already does not exist
    std::string directory = "output/" + out_filename + "/" + "Exact_Count";
    DIR* dir = opendir(directory.c_str());  // try to open the directory
    if (dir) {
        // Directory exists, go on to create a file in this location
        closedir(dir);
    }
    else if (ENOENT == errno){
        // Directory does not exist, try creating one
        std::string mkdir = "mkdir -p " + directory;
        if (std::system(mkdir.c_str()) == -1) {
            printf("Could not create directory output/%s/\n", out_filename.c_str());
            return;
        }

    }
    else {
        printf("Could not create/find directory output/%s/\n",out_filename.c_str());
        return;
    }

    out_filename = "output/" + out_filename + "/" + "Exact_Count" + "/" +
                   current_time + "-" + out_filename +   ".txt";
    FILE *f = fopen(out_filename.c_str(), "w");
    if (!f) {
        printf("Could not write output. Please check for write permissions. Lcaoltion: %s\n",out_filename.c_str());
        return;
    }

    fprintf(f, "#Filename = %s  \n", filename.c_str());
    fprintf(f, "########################\n");
    fprintf(f, "#Graph Properties\n");
    fprintf(f, "########################\n");
    fprintf(f, "vertices,edges,triangles\n");
    fprintf(f, "%lld,%lld,%lld\n", cg->nVertices, cg->nEdges, triangle_count);
    fclose(f);
}

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
    write_output(argv[1],&cg,out.estimate);
    return 0;
}