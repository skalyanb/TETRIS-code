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
    df_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1,2])
    median_error = df_error.iloc[0,0]
    max_error = df_error.iloc[0,1]

    return median_error,max_error


def plot_bar_median_max(file_name, data_dir, out_filename):
    Median_err = []
    Max_err = []

    for i,dir in enumerate(data_dir):
        for i, file in enumerate(os.listdir(dir)):
            data_file = dir + file
            if valid_file(dir,file):
                median_error, variance_error  = parse_outfile(data_file)
        Median_err.append(median_error)
        Max_err.append(variance_error)

    fig, ax = plt.subplots()

    barWidth = 0.2
    starting_pos = 1
    gap = 1
    minor_gap = 0.05

    med_pos = []
    max_pos = []
    # Set bar position
    for i in range(4):
        med_pos.append(starting_pos)
        max_pos.append(starting_pos+barWidth+minor_gap)
        starting_pos = starting_pos +1

    ax.bar(med_pos, Median_err, width=barWidth, color='b',  edgecolor= "black",label="Median Error %")
    ax.bar(max_pos, Max_err, width=barWidth, color='g',  edgecolor= "black",label="Max Error %")

    ax.legend(loc='upper right',fontsize=16)

    ax.axhline(y=2, color='r', linestyle='--')
    ax.axhline(y=5, color='r', linestyle='--')
    ax.text(1.02, 0.2, '2%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, 0.5, '5%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )

    y_min = 0
    y_max = 10
    ax.set_ylim(y_min,y_max)
    ax.tick_params(axis="y", labelsize=18)

    fig.suptitle('Accuracy of TETRIS',fontsize=20)
    ax.set_ylabel('Relative Error Percentage (%)',fontsize=20)
    ax.set_xlabel('Real world graph datasets',fontsize=20)

    ind = [1+i+0.2 for i in range(4)]
    algo_label = ['orkut','sinaweibo','twitter','friendster']
    ax.set_xticks(ind)
    ax.set_xticklabels(algo_label,fontsize=20)

    plt.rcParams.update(params)
    # plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # fig.set_size_inches(4,6)
    fig.savefig("output/plots/accuracy/"+out_filename+timestr+".eps",format='eps',bbox_inches="tight")



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
    # # title_info.append(f_name + ": 24M edges, 3M vertices")

    f_name = "soc-orkut"
    file_name.append(f_name)
    title_info.append(f_name + ": 213M edges")

    f_name = "soc-sinaweibo"
    file_name.append(f_name)
    title_info.append(f_name + ": 523M edges")
    #
    f_name = "soc-twitter-konect"
    file_name.append(f_name)
    title_info.append(f_name + ": 2.4B edges")
    #
    f_name = "soc-friendster"
    file_name.append(f_name)
    title_info.append(f_name + ": 3.6B edges")
    data_dir = []

    for i,file in enumerate(file_name):
        out_filename = "-accuracy-bar"
        data_dir.append("output/plot_data/intro_plot_data/"+file+".edges/")
    plot_bar_median_max(file_name, data_dir, out_filename)