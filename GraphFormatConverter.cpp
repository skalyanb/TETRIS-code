//
// Created by Suman Kalyan Bera on 2020-01-24.
//
#include <cstdio>
#include <chrono>

#include "include/Graph.h"
#include "include/GraphIO.h"

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cout << "Usage: ./GraphFormatConverter.out path/to/out/directory/ input_file" << std::endl;
        return 0;
    }

    Escape::Graph g;

    auto startTime = std::chrono::high_resolution_clock::now();
    if (loadGraph(argv[2], g, 1, Escape::IOFormat::escape))
        exit(1);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Loaded graph from "<< argv[2] <<"\nTime to load graph:" <<
              std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    Escape::CGraph cg = makeCSR(g);
    cg.sortById();
    endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Converted to CSR\n. Time taken:" <<
              std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << std::endl;

    // extract the input filename and construct output file path
    std::string in_filepath(argv[2]);
    std::string out_filename = in_filepath.substr(in_filepath.find_last_of("/\\") + 1);
    std::cout << out_filename << std::endl;
    std::string out_path(argv[1]);
    out_path = out_path + out_filename + ".csr";
    std::cout << out_path << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    cg.writeBinaryFile(out_path.c_str());
    endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Written to CSR format at" <<  out_path <<"\nTime taken:" <<
              std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << std::endl;

    return 0;
}