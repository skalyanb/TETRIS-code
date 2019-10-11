import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt
import statistics


dirname = os.path.dirname
sys.path.append(dirname(dirname(os.path.realpath(__file__))))

def skip_rows_err (index):
    if index <10 or index > 10:
        return True
    return False


def skip_rows_seen (index):
    if index <13 or index > 13:
        return True
    return False

def skip_rows (index):
    if index <4 or index > 5:
        return True
    return False

# Read edge_estimator file, read all the edge estimates and compute the median edge estimate.
# This will be written to a file, with the corresponding edge observe percentage
def read_edge_estimator (data_dir):
    edge_est = []
    edge_frac = 0.0
    actual_edge = 0.0
    for i, file in enumerate(os.listdir(data_dir)):
        path = os.path.join(data_dir, file)
        if os.path.isdir(path):
            # skip directories
            continue
        if not file.startswith('.'):
            data_file = data_dir + file
            df_tri = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x))
            actual_edge = df_tri["edges"][0]
            df_rawdata = pd.read_csv(data_file, sep=",",skiprows=15)
            edge_est.extend(df_rawdata["triangle_estimate"])
            df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
            edge_frac = df_edges_seen.iloc[0,0]
    return statistics.median(edge_est),edge_frac,actual_edge

def get_triangle_estimates(tri_cnt_dir):
    edges_seen=[]
    tri_est=[]
    triangle_count = 0
    for i, file in enumerate(os.listdir(tri_cnt_dir)):
        path = os.path.join(tri_cnt_dir, file)
        if os.path.isdir(path):
            # skip directories
            continue
        if not file.startswith('.'):
            data_file = tri_cnt_dir + file
            df_tri = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x))
            triangle_count = df_tri["triangles"][0]
            df_rawdata = pd.read_csv(data_file, sep=",",skiprows=15)
            tri_est.extend(df_rawdata["triangle_estimate"])
            edges_seen.extend(df_rawdata["fraction_of_edges_seen"])
    return edges_seen,tri_est, triangle_count

def edge_correction(edge_estimate_dir,tri_cnt_dir, out_dir, out_filename):
    edges_seen=[]
    tri_est=[]
    edge_est,edge_frac_offset,actual_edge = read_edge_estimator(edge_estimate_dir)
    tri_correction_offset = 1.0 * edge_est /  actual_edge

    for i, file in enumerate(os.listdir(tri_cnt_dir)):
        path = os.path.join(tri_cnt_dir, file)
        if os.path.isdir(path):
            # skip directories
            continue
        if not file.startswith('.'):
            data_file = tri_cnt_dir + file
            df_tri = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x))
            triangle_count = df_tri["triangles"][0]
            df_rawdata = pd.read_csv(data_file, sep=",",skiprows=15)
            tri_est = df_rawdata["triangle_estimate"]
            edges_seen = df_rawdata["fraction_of_edges_seen"]

            # Update the error statistics
            new_edges_seen = [x+edge_frac_offset for x in edges_seen]
            new_tri_est = [x*tri_correction_offset for x in tri_est]
            abs_err = [abs(x - triangle_count)*100.0/triangle_count for x in new_tri_est]
            stats = [statistics.mean(abs_err), statistics.median(abs_err), max(abs_err), statistics.stdev(abs_err)]

            # Write them to a new file.
            file = '_corrected_'+file
            out_path = os.path.join(out_dir, file)
            df_input_header = pd.read_csv(data_file, sep=",",nrows=16,header=None,names=[str(x) for x in range(10)])

            reader = open(path,'r')
            writer = open(out_path,'w')


            # fprintf(f, "#Filename = %s , Algo name = %s \n", params.filename.c_str(), params.algo_name.c_str());
            # fprintf(f, "########################\n");
            # fprintf(f, "#Graph Properties\n");
            # fprintf(f, "########################\n");
            # fprintf(f, "vertices,edges,triangles\n");
            # fprintf(f, "%lld,%lld,%lld\n", cg->nVertices, cg->nEdges, triangle_count);
            # fprintf(f, "no_of_repeat,seed_count,seed_vertex,degree_of_seed_vertex,walk_length,subsample_size,sparsification_prob\n");
            # fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lf\n",params.no_of_repeat, params.seed_count, params.seed_vertices[0],
                    # cg->degree(params.seed_vertices[0]),params.walk_length,
                    #     params.subsample_size, params.sparsification_prob);
            # fprintf(f, "#%s\n", algo_name.c_str());
            # fprintf(f, "#Results: Mean Err, Median Err, Max Err, stddev Err (in %%) of simple sampling\n");
            for i in range(0,10):
                line = reader.readline()
                writer.writelines(line)
                print (line)
            # fprintf(f, "%.3lf,%.3lf,%.3lf,%.3lf \n\n", est_stats.mean_error_percentage,
            #         est_stats.median_error_percentage, est_stats.max_error_percentage,
            #         est_stats.stddev_error_percentage);
            line = reader.readline()
            line = reader.readline()
            writer.write('%.3lf,%.3lf,%.3lf,%.3lf \n\n' % (stats[0],stats[1],stats[2],stats[3] ))
            # fprintf(f, "Fraction of edges seen, fraction of vertices seen(maximum over all run)\n");
            line = reader.readline()
            writer.writelines(line)
            # fprintf(f, "%.6lf,%.6lf\n\n", est_stats.edges_seen_max_percentage,
            # est_stats.vertices_seen_max_percentage);
            line = reader.readline()
            line = reader.readline()
            writer.write('%.6lf,0.0\n\n' % max(new_edges_seen))
            # fprintf(f, "triangle_estimate,fraction_of_edges_seen,fraction_of_vertices_seen\n");
            line = reader.readline()
            writer.writelines(line)
            # for (auto & est : estimates) {
            #     fprintf(f, "%.3lf,%.6lf,%.6lf\n", est.triangle_estimate, est.fraction_of_edges_seen,
            #             est.fraction_of_vertices_seen);
            # }
            for i,item in enumerate(new_tri_est):
                writer.write('%.3lf,%.6lf,0.0\n' % (new_tri_est[i],new_edges_seen[i]))
            reader.close()
            writer.close()

# There are three directory: 1. Edge estimate diretory which contain one single file with edge estimations
# 2. The directory with triangle counts that are not edge corrected.
# 3. The output directory where the new files will be written

if __name__ == "__main__":

    x_max = 1.1

    file_name = []
    title_info = []
    med_err = [] # Median error in edge count
    offset =[] # Percentage edges seen
    true_edge = [] # True edge count
    # f_name = "soc-flickr-und"
    # file_name.append(f_name)

    # title_info.append(f_name + ": 16M edges, 1.7M vertices")

    # f_name = "socfb-A-anon"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 24M edges, 3M vertices")
    #
    # f_name = "soc-orkut"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 106M edges, 3M vertices")
    #
    # f_name = "soc-sinaweibo"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges, 58M vertices")

    f_name = "soc-twitter-konect"
    file_name.append(f_name)
    title_info.append(f_name + ": 1.2B edges, 41M vertices")
    #
    # f_name = "soc-friendster"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")
    # #
    for i,file in enumerate(file_name):
        edge_estimate_dir = "output/plot_data/degree_bin_data/"+file+".edges/EstTriByRWandWghtedSampling/correction/"
        tri_cnt_dir = "output/plot_data/degree_bin_data/"+file+".edges/EstTriByRWandWghtedSampling/bin_3/"
        out_dir = "output/plot_data/degree_bin_data/"+file+".edges/EstTriByRWandWghtedSampling/edge_corrected/bin_3"
        out_filename = file

        edge_correction(edge_estimate_dir,tri_cnt_dir, out_dir,out_filename)

    #plot_comparison(data_dir_1,data_dir_2)

