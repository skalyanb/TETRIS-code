//
// Created by Suman Kalyan Bera on 2019-09-17.
// This file contains utility functions that analyze
// the performance of triangle estimators
//

#ifndef SUBGRAPHCOUNT_ESTIMATORUTIL_H
#define SUBGRAPHCOUNT_ESTIMATORUTIL_H

#include "TriangleEstimators.h"
#include "BaselineEstimators.h"

// TODO consider merging the below two structure into one
// A structure to store various statistics of an estimator across multiple runs

struct EstimatorStats {
    int no_of_repeats;
    double mean_error_percentage;
    double median_error_percentage;
    double stddev_error_percentage;
    double max_error_percentage;
    double vertices_seen_max_percentage;
    double edges_seen_max_percentage;
};


EstimatorStats GetErrorStatistics(std::vector<Estimates> algo_estimates, Count true_triangle_count) {

    int no_of_repeats = (int) algo_estimates.size();

    // Populate the error percentage for each run
    std::vector<double> algo_error_percentage_list(no_of_repeats);
    // TODO consider using lambda function
    for (int i = 0; i < no_of_repeats; i++)
        algo_error_percentage_list[i] =
                std::abs(algo_estimates[i].triangle_estimate - true_triangle_count) * 100 / true_triangle_count;

    int middle_index = no_of_repeats / 2;
    std::nth_element(algo_error_percentage_list.begin(), algo_error_percentage_list.begin() + middle_index,
                     algo_error_percentage_list.end());
    double median_error_percentage = algo_error_percentage_list[middle_index];

    double max_error_percentage = *std::max_element(algo_error_percentage_list.begin(),
                                                    algo_error_percentage_list.end());

    double sum = std::accumulate(algo_error_percentage_list.begin(), algo_error_percentage_list.end(), 0.0);
    double mean_error_percentage = sum / no_of_repeats;

    std::vector<double> diff(no_of_repeats);
    std::transform(algo_error_percentage_list.begin(), algo_error_percentage_list.end(), diff.begin(),
                   [mean_error_percentage](double x) { return x - mean_error_percentage; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev_error_percentage = std::sqrt(sq_sum / no_of_repeats);


    // Maximum edges and vertices seen over all the runs
    double vertices_seen_percentage = (*std::max_element(algo_estimates.begin(),
                                                         algo_estimates.end(),
                                                         ComparatorByVerticesSeen)).fraction_of_vertices_seen;
    double edges_seen_percentage = (*std::max_element(algo_estimates.begin(),
                                                      algo_estimates.end(),
                                                      ComparatorByEdgesSeen)).fraction_of_edges_seen;

    EstimatorStats est_stats = {no_of_repeats, mean_error_percentage, median_error_percentage,
                              stdev_error_percentage, max_error_percentage,
                              vertices_seen_percentage, edges_seen_percentage};
    return est_stats;
}


// TODO use c++ file buffer here
void WriteHeaderInOutput(FILE *f, Parameters params, CGraph *cg, Count triangle_count) {
    fprintf(f, "#Filename = %s , Algo name = %s \n", params.filename.c_str(), params.algo_name.c_str());
    fprintf(f, "########################\n");
    fprintf(f, "#Graph Properties\n");
    fprintf(f, "########################\n");
    fprintf(f, "vertices,edges,triangles\n");
    fprintf(f, "%lld,%lld,%lld\n", cg->nVertices, cg->nEdges, triangle_count);
    fprintf(f, "no_of_repeat,seed_count,walk_length,subsample_size,sparsification_prob\n");
    fprintf(f, "%d,%lld,%lld,%lld,%.3lf\n",params.no_of_repeat, params.seed_count, params.walk_length,
            params.subsample_size, params.sparsification_prob);

}

void WriteAlgorithmOutput(FILE *f, std::string algo_name, Parameters params,
                          EstimatorStats est_stats) {
    fprintf(f, "#%s\n", algo_name.c_str());
    fprintf(f, "#Results: Mean Err, Median Err, Max Err, stddev Err (in %%) of simple sampling\n");
    fprintf(f, "%.3lf,%.3lf,%.3lf,%.3lf \n\n", est_stats.mean_error_percentage,
            est_stats.median_error_percentage, est_stats.max_error_percentage,
            est_stats.stddev_error_percentage);
    fprintf(f, "Fraction of edges seen, fraction of vertices seen(maximum over all run)\n");
    fprintf(f, "%.3lf,%.3lf\n\n", est_stats.edges_seen_max_percentage,
            est_stats.vertices_seen_max_percentage);
}

void WriteRawData (FILE *f, std::vector<Estimates> const &estimates) {
    fprintf(f, "triangle_estimate,fraction_of_edges_seen,fraction_of_vertices_seen\n");
    for (auto & est : estimates) {
        fprintf(f, "%.3lf,%.3lf,%.3lf\n", est.triangle_estimate, est.fraction_of_edges_seen,
                est.fraction_of_vertices_seen);
    }
}

std::string GetTimestamp() {
    auto now = std::time(nullptr);
    char buf[sizeof("YYYY-MM-DD_HH:MM:SS")];
    return std::string(buf, buf +
                            std::strftime(buf, sizeof(buf), "%F_%T", std::gmtime(&now)));
}


template <typename TF>
void TriangleEstimator (CGraph *cg, Parameters params, Count true_triangle_count, TF func) {
    std::vector<Estimates> estimates;

    for (Count i = 0; i < params.no_of_repeat; i++) {
        Estimates est = func(cg, params);
        estimates.push_back(est);
        std::cout << i << "\n";
    }
    EstimatorStats est_stats = GetErrorStatistics(estimates, true_triangle_count);

    // Print to console
//    WriteHeaderInOutput(stdout, params, cg, true_triangle_count);
//    WriteAlgorithmOutput(stdout, params.algo_name, params, est_stats);
//    WriteRawData(stdout,estimates);


    // print to file
    std::string current_time = GetTimestamp();
    std::string out_filename = params.filename.substr(params.filename.find_last_of("/\\") + 1);
    current_time = "output/" + out_filename +"/" +
                current_time + "-" + out_filename +  "-" + params.algo_name + "-" +
                std::to_string(params.sparsification_prob) +".txt";
    FILE *f = fopen(current_time.c_str(), "w");
    if (!f) {
        printf("could not write to output to out.txt\n");
        return;
    }
    WriteHeaderInOutput(f, params, cg, true_triangle_count);
    WriteAlgorithmOutput(f, params.algo_name, params, est_stats);
    WriteRawData(f,estimates);


    fclose(f);


}

#endif //SUBGRAPHCOUNT_ESTIMATORUTIL_H
