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
    edges = []
    medians=[]
    vertices=[]
    query = []
    for i, file in enumerate(os.listdir(data_dir)):
        data_file = data_dir + file
        if valid_file(data_dir,file):
            median_error, edges_seen, vertices_seen, query_made = parse_outfile(data_file)
            medians.append(median_error)
            edges.append(edges_seen)
            vertices.append(vertices_seen)
            query.append(query_made)

    return edges,medians,vertices,query

def make_plot_edge (X,Y,filename, title_info, algo_list, marker_list, color_list, save):

    #scatter plot
    fig, ax_edge = plt.subplots()

    for i,algo in enumerate(algo_list):
        ax_edge.plot(X[i], Y[i], c=color_list[i], marker=marker_list[i], label=algo, markersize=12)
    # ax_edge.plot(X[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)
    # ax_edge.plot(X[1], Y[1], c='red', marker='o', label='SRW1',markersize=12)
    # # ax_edge.plot(X[2], Y[2], c='red', marker='^', label='SERWC',markersize=12)
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
    #plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps',bbox_inches="tight")



def make_plot_query (X,Y,filename, title_info, algo_list, marker_list, color_list, save):

    #scatter plot
    fig, ax_edge = plt.subplots()

    for i,algo in enumerate(algo_list):
        if algo == "TETRIS":
            ax_edge.plot(X[i], Y[i], c=color_list[i], marker=marker_list[i], label=algo, linewidth=3, markersize=12)
        else:
            ax_edge.plot(X[i], Y[i], c=color_list[i], marker=marker_list[i], label=algo, markersize=12)    # ax_edge.plot(X[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)

    ax_edge.set_xlim(0.15,1.7)
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
    #plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/baseline_comparison/query/basu/"+filename+timestr+".eps",format='eps',bbox_inches="tight")


def make_plot_vertex(Z,Y,filename, title_info, algo_list, marker_list, color_list, save):

    #scatter plot
    fig, ax_vertex = plt.subplots()

    for i,algo in enumerate(algo_list):
        if algo == "TETRIS":
            ax_vertex.plot(Z[i], Y[i], c=color_list[i], marker=marker_list[i], label=algo, linewidth=10, markersize=12)
        else:
            ax_vertex.plot(Z[i], Y[i], c=color_list[i], marker=marker_list[i], label=algo, markersize=12)
    # ax_vertex.plot(Z[1], Y[1], c='red', marker='o', label='SEC',markersize=12)
    # ax_vertex.plot(Z[2], Y[2], c='red', marker='^', label='SERWC',markersize=12)
    # ax_vertex.plot(Z[3], Y[3], c='red', marker='D', label='UESS',markersize=12)
    # ax_vertex.plot(Z[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)
    # ax_vertex.set_xlim(1,8)
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
    #plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps',bbox_inches="tight")

def plot_comparison(data_dir, filename, title_info, algo_list, marker_list, color_list, save=True):

    edges_list = []
    medians_list= []
    vertices_list = []
    query_list = []
    for i,dir in enumerate(data_dir):
        tempEdges, tempMedians, tempVertices, tempQuery = populate_data(dir)
        edges_list.append(tempEdges)
        medians_list.append(tempMedians)
        vertices_list.append(tempVertices)
        query_list.append(tempQuery)
    # make_plot_edge (edges_list,medians_list,filename, title_info, algo_list, marker_list, color_list, save)
    # make_plot_vertex(vertices_list,medians_list,filename, title_info, algo_list, marker_list, color_list, save)

def make_plot_query_combined (query,medians,filename, title_info, algo_list, marker_list, color_list, save):
    #scatter plot
    fig, ax_edge = plt.subplots(1,4, sharex=True)

    plt.subplots_adjust(bottom = 0.19, top = 0.8, right = 0.95, left = 0.1, hspace = 0.38)

    for ind, file in enumerate(filename):
        ax_edge[ind].set_title(title_info[ind],fontsize=14)
        for i,algo in enumerate(algo_list):
            if algo == "TETRIS":
                ax_edge[ind].plot(query[file][i], medians[file][i], c=color_list[i], marker=marker_list[i], label=algo, linewidth=3, markersize=8)
            else:
                ax_edge[ind].plot(query[file][i], medians[file][i], c=color_list[i], marker=marker_list[i], label=algo, markersize=8)    # ax_edge.plot(X[0], Y[0], c='blue', marker='s', label='TETRIS',markersize=12)

    # ax_edge[1].set_xlabel('Percentage of queries %',fontsize=12)
    fig.text(0.52,0.11,"Percentage of queries %",ha="center",va="center",fontsize=14)
    ax_edge[0].set_ylabel('Median Relative Error %' ,fontsize=14)

    for ax in fig.get_axes():
        ax.set_xlim(0.15,1.75)
        ax.set_ylim(0,10)
        # ax_edge.legend(loc='upper right',fontsize=18)
        #add x and y labels

        ax.tick_params(axis="y", labelsize=12)
        ax.tick_params(axis="x", labelsize=12)

        # ax.set_aspect(0.2)
        # ax.set_adjustable('box')

        # ax.label_outer()

    # Set figure title and axis labels
    handles, labels = ax_edge[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol = 5, fontsize = 12)

#add title
    #fig.suptitle(title_info,fontsize=18)
    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if (save):
        fig.savefig("output/plots/baseline_comparison/query/basu/combined"+timestr+".eps",format='eps',bbox_inches="tight")

def plot_comparison_combined(input_dir, file_name, title_info, algo_list, marker_list, color_list, save=True):

    edges = {}
    medians= {}
    vertices = {}
    query = {}

    for i,file in enumerate(file_name):
        edges_list = []
        medians_list= []
        vertices_list = []
        query_list = []

        for j,algo in enumerate(algo_list):
            data_dir = input_dir+file+".edges.csr/"+algo+"/"
            tempEdges, tempMedians, tempVertices, tempQuery = populate_data(data_dir)
            edges_list.append(tempEdges)
            medians_list.append(tempMedians)
            vertices_list.append(tempVertices)
            query_list.append(tempQuery)

        edges [file] = edges_list
        medians [file] = medians_list
        vertices [file] = vertices_list
        query [file] = query_list
    # make_plot_edge (edges_list,medians_list,filename, title_info, algo_list, marker_list, color_list, save)
    # make_plot_vertex(vertices_list,medians_list,filename, title_info, algo_list, marker_list, color_list, save)
    make_plot_query_combined (query,medians,file_name, title_info, algo_list, marker_list, color_list, save=True)

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
    #title_info.append(f_name + ": 31M edges, 3M vertices")
    title_info.append("orkut")

    # # title_info.append(f_name + ": 3M vertices")

    f_name = "soc-sinaweibo"
    file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges, 59M vertices")
    title_info.append("weibo")

    f_name = "soc-twitter-konect"
    file_name.append(f_name)
    # title_info.append(f_name + ": 1.2B edges, 41M vertices")
    title_info.append("twitter")


    f_name = "soc-friendster"
    file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")
    title_info.append("friendster")

    # Potential list of algorithms: TETRIS,SRW1,VertexMCMC,UESS,SEC,SERWC
    # algo_list = ["TETRIS","SRW1","VertexMCMC","UESS","SEC","SERWC"]
    algo_list = ["TETRIS","SRW1","VertexMCMC", "SERWC", "RWS"]
    marker_list = ['s','o','D','*','^','p']
    color_list = ['blue','green','magenta','red','red','red']
    save = True
    input_dir = "output/plot_data/baseline_TETRIS_SRW_MCMC/from_basu/"

    # for i,file in enumerate(file_name):
    #     data_dir = []
    #     out_filename = file + "-comparison-"
    #     for j,algo in enumerate(algo_list):
    #         data_dir.append(input_dir+file+".edges.csr/"+algo+tag)
    #     plot_comparison(data_dir, out_filename, title_info[i], algo_list, marker_list, color_list, save)

    plot_comparison_combined(input_dir, file_name, title_info, algo_list, marker_list, color_list, save)