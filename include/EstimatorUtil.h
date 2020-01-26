//
// Created by Suman Kalyan Bera on 2019-09-17.
// This file contains utility functions that analyze
// the performance of triangle estimators
//

#ifndef SUBGRAPHCOUNT_ESTIMATORUTIL_H
#define SUBGRAPHCOUNT_ESTIMATORUTIL_H

#include <cstdlib>
#include <chrono>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>


#include "EstimatorUtilStruct.h"
#include "EstimatorUtilStats.h"
#include "TriangleEstimators.h"
//#include "BaselineEstimators.h"
#include "baseline/VertexMCMC.h"
#include "baseline/SubgraohRandomWalk_SRW.h"
#include "TETRIS.h"


template <typename TF>
void TriangleEstimator (CGraph *cg, Parameters params, Count true_triangle_count, TF func) {
    std::vector<Estimates> estimates;

    for (Count i = 0; i < params.no_of_repeat; i++) {
        auto startTime = std::chrono::high_resolution_clock::now();
        Estimates est = func(cg, params);
        auto endTime = std::chrono::high_resolution_clock::now();
        std::cout << "Run " << i<<". Time taken for this run = " <<
                  std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()
                  << " seconds" << std::endl;
        estimates.push_back(est);
    }
    EstimatorStats est_stats = GetErrorStatistics(estimates, true_triangle_count);

    // Print to console
    if (params.print_to_console) {
        WriteHeaderInOutput(stdout, params, cg, true_triangle_count);
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
        WriteHeaderInOutput(f, params, cg, true_triangle_count);
        WriteAlgorithmOutput(f, params.algo_name, params, est_stats);
        WriteRawData(f, estimates);

        fclose(f);
    }

}

void EdgeEstimatorUtil (CGraph *cg, Parameters params, int c) {
    OrderedEdgeCollection randomEdgeCollection;
    std::vector<Estimates> estimates;
    // Set up random number generator
    std::random_device rd;

    for (Count i = 0; i < params.no_of_repeat; i++) {
        // Using this random number generator initializize a PRNG: this PRNG is passed along to
        // draw an element from various distribution
        std::mt19937 mt(rd());

        randomEdgeCollection = GetEdgesByRandomWalk(cg, params, mt);
        Estimates est = EstimateEdgeCount(cg, randomEdgeCollection, params, c);
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
                       std::to_string(params.sparsification_prob) + "-" + std::to_string(c) +
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

void EdgeEstimator (CGraph *cg, Parameters params) {
    std::vector<int> c = {10};
    for (auto item : c) {
        EdgeEstimatorUtil(cg, params, item);
    }
}


#endif //SUBGRAPHCOUNT_ESTIMATORUTIL_H
