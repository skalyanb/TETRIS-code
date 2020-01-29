//
// Created by Suman Kalyan Bera on 2019-10-09.
//

#ifndef SUBGRAPHCOUNT_ESTIMATORUTILSTATS_H
#define SUBGRAPHCOUNT_ESTIMATORUTILSTATS_H

#include <cstdlib>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <ctime>

#include "EstimatorUtilStruct.h"

EstimatorStats GetErrorStatistics(std::vector<Estimates> algo_estimates, Count true_count) {

    int no_of_repeats = (int) algo_estimates.size();

    // Populate the error percentage for each run
    std::vector<double> algo_error_percentage_list(no_of_repeats);
    // TODO consider using lambda function
    for (int i = 0; i < no_of_repeats; i++)
        algo_error_percentage_list[i] =
                std::abs(algo_estimates[i].estimate - true_count) * 100 / true_count;

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

    double query_complexity_percentage = (*std::max_element(algo_estimates.begin(),
                                                      algo_estimates.end(),
                                                      ComparatorByEdgesSeen)).query_complexity;

    EstimatorStats est_stats = {no_of_repeats, mean_error_percentage, median_error_percentage,
                                stdev_error_percentage, max_error_percentage,
                                vertices_seen_percentage, edges_seen_percentage, query_complexity_percentage};
    return est_stats;
}


// TODO use c++ file buffer here
void WriteHeaderInOutput(FILE *f, Parameters params, CGraph *cg, Count triangle_count) {
    fprintf(f, "#Filename = %s , Algo name = %s, "
               "edge count available=%s, CSS = %s, NB=%s \n",
            params.filename.c_str(), params.algo_name.c_str(),
            params.normalization_count_available ? "true":"false",
            params.CSS ? "true":"false",
            params.NB ? "true":"false");
    fprintf(f, "########################\n");
    fprintf(f, "#Algorithm : %s\n",params.algo_name.c_str());
    fprintf(f, "########################\n");
    fprintf(f, "vertices,edges,triangles\n");
    fprintf(f, "%lld,%lld,%lld\n", cg->nVertices, cg->nEdges, triangle_count);
    fprintf(f, "no_of_repeat,seed_count,seed_vertex,degree_of_seed_vertex,"
               "walk_length,subsample_size,sparsification_prob\n");
    fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lf\n",params.no_of_repeat, params.seed_count, params.seed_vertices[0],
            cg->degree(params.seed_vertices[0]),params.walk_length,
            params.subsample_size, params.sparsification_prob);

}

void WriteAlgorithmOutput(FILE *f, std::string algo_name, Parameters params,
                          EstimatorStats est_stats) {
    fprintf(f, "#--------------------------------------------------------%s\n", algo_name.c_str());
    fprintf(f, "#Results: Mean Err, Median Err, Max Err, stddev Err (in %%) of simple sampling\n");
    fprintf(f, "%.3lf,%.3lf,%.3lf,%.3lf \n\n", est_stats.mean_error_percentage,
            est_stats.median_error_percentage, est_stats.max_error_percentage,
            est_stats.stddev_error_percentage);
    fprintf(f, "Query Complexity, Fraction of edges seen, fraction of vertices seen(maximum over all run)\n");
    fprintf(f, "%.6lf,%.6lf,%.6lf\n\n", est_stats.query_complexity_max_percentage, est_stats.edges_seen_max_percentage,
            est_stats.vertices_seen_max_percentage);
}

void WriteRawData (FILE *f, std::vector<Estimates> const &estimates) {
    fprintf(f, "triangle_estimate,fraction_of_edges_seen,fraction_of_vertices_seen\n");
    for (auto & est : estimates) {
        fprintf(f, "%.3lf,%.6lf,%.6lf\n", est.estimate, est.fraction_of_edges_seen,
                est.fraction_of_vertices_seen);
    }
}

std::string GetTimestamp() {
    auto now = std::time(nullptr);
    char buf[sizeof("YYYY-MM-DD_HH:MM:SS")];
    return std::string(buf, buf +
                            std::strftime(buf, sizeof(buf), "%F_%T", std::gmtime(&now)));
}



#endif //SUBGRAPHCOUNT_ESTIMATORUTILSTATS_H
