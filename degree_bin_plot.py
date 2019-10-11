import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt
import statistics
import numpy as np



dirname = os.path.dirname
sys.path.append(dirname(dirname(os.path.realpath(__file__))))

def skip_rows_err (index):
    if index <10 or index > 10:
        return True
    return False

def plot_degree_bin (bin_dir, out_filename, title_info):
    # We will go into each directory in the bin_dir list
    # and then iterate over each file inside the directory.
    # For each file we take the median accuracy and the variance accuracy
    #
    median_err=[]
    variance_err=[]


    # Enumerate over each bin
    for i,data_dir in enumerate(bin_dir):
        # We are inside a bin; Go over each file and collect the median accuracy and variance.
        bin_median_error = []
        bin_variance_error = []
        for i, file in enumerate(os.listdir(data_dir)):
            path = os.path.join(data_dir, file)
            if os.path.isdir(path):
                # skip directories
                continue
            if not file.startswith('.'):
                data_file = data_dir + file
                df_median_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1,3])
                bin_median_error.append(df_median_error.iloc[0,0])
                bin_variance_error.append(df_median_error.iloc[0,1])
        median_err.append(bin_median_error)
        variance_err.append(bin_variance_error)

    num_bin = len(bin_dir)
    group_len = len(median_err[0])

    # Set bar width
    barWidth = 0.2
    starting_pos = 1
    gap = 1

    x={}
    # Set bar position
    for i in range(num_bin):
        x[i] = [starting_pos+i*barWidth for i in range(group_len)]
        starting_pos = starting_pos + gap + barWidth*group_len



    fig, ax = plt.subplots()
    ax.bar(x[0], median_err[0], width=barWidth, color='b',  edgecolor= "black",label="Gender")
    ax.bar(x[1], median_err[1], width=barWidth, color='b',  edgecolor= "black",label="Type")
    ax.bar(x[2], median_err[2], width=barWidth, color='b',  edgecolor= "black",label="Type")
    ax.bar(x[3], median_err[3], width=barWidth, color='b',  edgecolor= "black",label="Type")
    #rects4 = ax.bar(x_4, median_err[4], width=barWidth, color='b',  edgecolor= "black",label="Type")

    y_min = 0
    y_max = 10
    ax.set_ylim(y_min,y_max)

    # y_ind = [i/2 for i in range(4)]
    # y_label = [str(i/2)+'%' for i in range(4)]
    # ax.set_yticks(y_ind)
    # ax.set_yticklabels(y_label)

    ax.set_title('Robustness of TETRIS\n'+title_info,fontsize=16)

    ax.set_ylabel('Absolute Median Error (%)',fontsize=16)

    ind = [x[i][1] for i in range(num_bin)]
    ax.set_xticks(ind)
    ax.set_xticklabels(('deg$\in [1,10]$', 'deg$\in [10,10^2]$','deg$\in [10^2,10^3]$', 'deg$\in [10^3,10^4]$'),fontsize=11)
    ax.set_xlabel('Degree based buckets',fontsize=16)
    plt.show()

    timestr = time.strftime("%Y%m%d-%H%M%S")

    fig.savefig("output/plots/degree-bin/"+out_filename+"-"+timestr+".eps",format='eps')




if __name__ == "__main__":

    x_max = 1.1

    file_name = []
    title_info = []

    # f_name = "soc-flickr-und"
    # file_name.append(f_name)

    # title_info.append(f_name + ": 16M edges, 1.7M vertices")

    f_name = "socfb-A-anon"
    file_name.append(f_name)
    title_info.append(f_name + ": 24M edges")
    no_bin = 4
    # #
    f_name = "soc-orkut"
    file_name.append(f_name)
    title_info.append(f_name + ": 106M edges")

    # f_name = "soc-sinaweibo"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges, 58M vertices")
    #
    f_name = "soc-twitter-konect"
    file_name.append(f_name)
    title_info.append(f_name + ": 1.2B edges")
    no_bin = 4

    # f_name = "soc-friendster"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")
    # #
    for i,file in enumerate(file_name):
        bin_directory = []
        # The main data dirextory
        data_dir = "output/plot_data/degree_bin_data/"+file+".edges/EstTriByRWandWghtedSampling/bin"
        for j in range(no_bin):
            bin_directory.append(data_dir+'_'+str(j)+'/')
        out_filename = file
        plot_degree_bin(bin_directory,out_filename,title_info[i])


