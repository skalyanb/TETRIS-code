//
// Created by Suman Kalyan Bera on 2019-10-08.
//

//
// Created by Suman Kalyan Bera on 2019-09-28.
//

#include <iostream>
#include <cstdlib>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>


#include "include/GraphIO.h"
#include "include/Graph.h"
//#include "include/TriangleEstimators.h"
//#include "include/EstimatorUtil.h"

using namespace Escape;

std::string GetTimestamp() {
    auto now = std::time(nullptr);
    char buf[sizeof("YYYY-MM-DD_HH:MM:SS")];
    return std::string(buf, buf +
                            std::strftime(buf, sizeof(buf), "%F_%T", std::gmtime(&now)));
}

void write_output (std::string filename, CGraph *cg, Count sum_deg) {

    // take current timestamp
    std::string current_time = GetTimestamp();
    // extract the input filename
    std::string out_filename = filename.substr(filename.find_last_of("/\\") + 1);

    // create the directory output/filename/algoname, if it already does not exist
    std::string directory = "output/" + out_filename + "/" + "Properties";
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

    out_filename = "output/" + out_filename + "/" + "Properties" + "/" +
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
    fprintf(f, "vertices,edges,sum_degree\n");
    fprintf(f, "%lld,%lld,%lld\n", cg->nVertices, cg->nEdges, sum_deg);
    fclose(f);
}

int main(int argc, char *argv[]) {

    if (argc!=2) {
        std::cout << "Usage: ./GraphProperties.out input_file_name";
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

    Count sum_degree =0;
    for (VertexIdx src; src < cg.nVertices; src++) {
        VertexIdx src_deg = cg.degree(src);
        for (EdgeIdx start = cg.offsets[src]; start < cg.offsets[src+1]; start++) {
            VertexIdx dst = cg.nbors[start];
            VertexIdx dst_deg = cg.degree(dst);
            // {src,dst} is the edge
            if ( src_deg < dst_deg)
                sum_degree += src_deg;
            else
                sum_degree += dst_deg;
        }
    }

    Count sum_degree_check = 0;
    VertexIdx n = cg.nVertices;
    for (EdgeIdx e=0; e < cg.nEdges; e++){
        VertexIdx dst = cg.nbors[e];
        VertexIdx src = std::upper_bound (cg.offsets, cg.offsets+n, e) - cg.offsets -1;
        VertexIdx src_deg = cg.degree(src);
        VertexIdx dst_deg = cg.degree(dst);
        // {src,dst} is the edge
        if ( src_deg < dst_deg)
            sum_degree_check += src_deg;
        else
            sum_degree_check += dst_deg;
    }
    if (sum_degree!=sum_degree_check) {
        printf("Mismatch! Error. Stop. Quit");
        exit(-1);
    }
    write_output(argv[1],&cg,sum_degree);
    return 0;
}