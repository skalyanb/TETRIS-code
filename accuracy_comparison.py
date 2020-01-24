import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt

params = {'legend.fontsize': 'xx-large',
          'axes.labelsize': 'xx-large',
          'axes.titlesize':'xx-large',
          'xtick.labelsize':'xx-large',
          'ytick.labelsize':'xx-large'}


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

def valid_file (data_dir,file):
    path = os.path.join(data_dir, file)
    if os.path.isdir(path):
        # skip directories
        return False
    if file.startswith('.'):
        return False
    return True


def parse_outfile (data_file):
    df_median_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1])
    median_error = df_median_error.iloc[0,0]

    df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
    edges_seen = df_edges_seen.iloc[0,0]

    # df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[1])
    # vertices_seen = df_edges_seen.iloc[0,0]

    return median_error,edges_seen

def populate_data (data_dir):
    median_err_list = []
    edges_seen_list=[]
    for i, file in enumerate(os.listdir(data_dir)):
        data_file = data_dir + file
        if valid_file(data_dir,file):
            median_error, edges_seen = parse_outfile(data_file)
            median_err_list.append(median_error)
            edges_seen_list.append(edges_seen)
            # Z.append(vertices_seen)

    return median_err_list,edges_seen_list


def plot_accuracy_varying(file_name, data_dir):

    median = []
    edges = []

    for i,dir in enumerate(data_dir):
        temp_median_err, temp_edges_seen = populate_data(dir)
        median.append(temp_median_err)
        edges.append(temp_edges_seen)

    #scatter plot
    fig, ax = plt.subplots()

    marker_list = ['s','o','v','s','D']
    color_list = ['b','g','r','y','m']
    for i,file in enumerate(file_name):
        ax.plot(edges[i], median[i], c=color_list[i], marker=marker_list[i], label=file)
    ax.set_xlim(2,3.3)
    ax.set_ylim(0,10)
    ax.tick_params(axis="y", labelsize=18)
    ax.tick_params(axis="x", labelsize=18)
    ax.legend(loc='upper right',fontsize=16)

    #add x and y labels
    ax.set_xlabel('Percentage of Edges Visited (%)',fontsize=20)
    ax.set_ylabel('Median Relative Error (%)',fontsize=20)
    fig.suptitle('Accuracy of TETRIS',fontsize=20)

    #show plot
    plt.rcParams.update(params)
    # plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # fig.set_size_inches(4,6)
    fig.savefig("output/plots/accuracy/varying_walk_length"+timestr+".eps",format='eps',bbox_inches="tight")



if __name__ == "__main__":


    file_name = []
    title_info = []

    tag ="/"
    #tag = "/vertex/"

    # f_name = "soc-flickr-und"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 16M edges, 1.7M vertices")

    # f_name = "socfb-A-anon"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 24M edges, 3M vertices")

    f_name = "soc-orkut"
    file_name.append(f_name)
    # title_info.append(f_name + ": 106M edges, 3M vertices")

    f_name = "soc-sinaweibo"
    file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges, 58M vertices")
    #
    f_name = "soc-twitter-konect"
    file_name.append(f_name)
    # title_info.append(f_name + ": 1.2B edges, 41M vertices")
    #
    f_name = "soc-friendster"
    file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")
    data_dir = []
    for i,file in enumerate(file_name):
        data_dir.append("output/plot_data/accuracy_plot_data/varying_walk_length/"+file+".edges/")

    plot_accuracy_varying(file_name, data_dir)