import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt

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

    df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[1])
    edges_seen = df_edges_seen.iloc[0,0]

    df_vertices_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[2])
    vertices_seen = df_vertices_seen.iloc[0,0]

    df_query = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
    query = df_query.iloc[0,0]

    return median_error,edges_seen, vertices_seen, query

def populate_data (data_dir):
    X = []
    Y=[]
    Z=[]
    W = []
    for i, file in enumerate(os.listdir(data_dir)):
        data_file = data_dir + file
        if valid_file(data_dir,file):
            median_error, edges_seen, vertices_seen, query = parse_outfile(data_file)
            Y.append(median_error)
            X.append(edges_seen)
            Z.append(vertices_seen)
            W.append(query)

    return X,Y,Z,W

def make_plot_edge (X,Y,filename, title_info, save):

    #scatter plot
    fig, ax_edge = plt.subplots()

    ax_edge.plot(X[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)
    ax_edge.plot(X[1], Y[1], c='red', marker='o', label='SRW1',markersize=12)
    # ax_edge.plot(X[2], Y[2], c='red', marker='^', label='SERWC',markersize=12)
    # ax_edge.plot(X[3], Y[3], c='red', marker='D', label='UESS',markersize=12)
    #ax_edge.set_xlim(0.15,1.1)
    ax_edge.set_ylim(0,20)
    ax_edge.legend(loc='upper right',fontsize=18)
    #add x and y labels
    ax_edge.set_xlabel('Percentage of Edges Visited',fontsize=20)
    ax_edge.set_ylabel('Median Relative Error %' ,fontsize=20)

    ax_edge.tick_params(axis="y", labelsize=18)
    ax_edge.tick_params(axis="x", labelsize=18)
    #add title
    fig.suptitle(title_info,fontsize=18)
    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps',bbox_inches="tight")



def make_plot_query (X,Y,filename, title_info, save):

    #scatter plot
    fig, ax_edge = plt.subplots()

    ax_edge.plot(X[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)
    ax_edge.plot(X[1], Y[1], c='red', marker='o', label='SRW1',markersize=12)
    # ax_edge.plot(X[2], Y[2], c='red', marker='^', label='SERWC',markersize=12)
    # ax_edge.plot(X[3], Y[3], c='red', marker='D', label='UESS',markersize=12)
    #ax_edge.set_xlim(0.15,1.1)
    ax_edge.set_ylim(0,20)
    ax_edge.legend(loc='upper right',fontsize=18)
    #add x and y labels
    ax_edge.set_xlabel('Percentage of queries %',fontsize=20)
    ax_edge.set_ylabel('Median Relative Error %' ,fontsize=20)

    ax_edge.tick_params(axis="y", labelsize=18)
    ax_edge.tick_params(axis="x", labelsize=18)
    #add title
    fig.suptitle(title_info,fontsize=18)
    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/baseline_comparison/"+filename+timestr+".eps",format='eps',bbox_inches="tight")


def make_plot_vertex(Z,Y,filename, title_info, save):

    #scatter plot
    fig, ax_vertex = plt.subplots()

    ax_vertex.plot(Z[1], Y[1], c='red', marker='o', label='SEC',markersize=12)
    ax_vertex.plot(Z[2], Y[2], c='red', marker='^', label='SERWC',markersize=12)
    ax_vertex.plot(Z[3], Y[3], c='red', marker='D', label='UESS',markersize=12)
    ax_vertex.plot(Z[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)
    ax_vertex.set_xlim(1,8)
    ax_vertex.set_ylim(0,15)
    ax_vertex.legend(loc='upper right',fontsize=18)
    #add x and y labels
    ax_vertex.set_xlabel('Percentage of Vertices Visited',fontsize=20)
    ax_vertex.set_ylabel('Median Relative Error %',fontsize=20)


    ax_vertex.tick_params(axis="y", labelsize=18)
    ax_vertex.tick_params(axis="x", labelsize=18)

    #add title
    fig.suptitle(title_info,fontsize=18)

    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps',bbox_inches="tight")


def plot_comparison(data_dir, filename, title_info, save=True):

    X = []
    Y= []
    Z = []
    W = []

    for i,dir in enumerate(data_dir):
        tempX, tempY, tempZ, tempW = populate_data(dir)
        X.append(tempX)
        Y.append(tempY)
        Z.append(tempZ)
        W.append(tempW)
    # make_plot_edge (X,Y,filename, title_info, save)
    # make_plot_vertex(Z,Y,filename, title_info, save)
    make_plot_query (W,Y,filename, title_info, save)

if __name__ == "__main__":


    file_name = []
    title_info = []

    tag ="/"
    # tag = "/vertex/"

    # f_name = "soc-flickr-und"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 16M edges, 1.7M vertices")

    # f_name = "socfb-A-anon"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 24M edges, 3M vertices")

    f_name = "soc-orkut"
    file_name.append(f_name)
    # title_info.append(f_name + ": 3M vertices")
    title_info.append(f_name + ": 31M edges")

    # f_name = "soc-sinaweibo"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges")
    #
    # f_name = "soc-twitter-konect"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.2B edges, 41M vertices")
    #
    # f_name = "soc-friendster"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")

    algo_list = []

    algo_list.append("TETRIS")
    algo_list.append("SRW1")

    for i,file in enumerate(file_name):
        data_dir = []
        out_filename = file + "-comparison-"
        for j,algo in enumerate(algo_list):
            data_dir.append("output/plot_data/new_baseline/"+file+".edges.csr/"+algo+tag)
        plot_comparison(data_dir, out_filename, title_info[i])