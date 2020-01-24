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
    df_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1,3])
    median_error = df_error.iloc[0,0]
    variace_error = df_error.iloc[0,1]

    return median_error,variace_error


def plot_accuracy(file_name, data_dir, out_filename, title_info=" "):

    Median = []
    Variance = []

    for i,dir in enumerate(data_dir):
        for i, file in enumerate(os.listdir(dir)):
            data_file = dir + file
            if valid_file(dir,file):
                median_error, variance_error  = parse_outfile(data_file)
        Median.append(median_error)
        Variance.append(variance_error)

    fig, ax = plt.subplots()

    x_axis = list(range(1,len(file_name)+1))
    algo_label = ['orkut','sinaweibo','twitter','friendster']
    ax.errorbar(x_axis, Median, yerr=Variance, fmt='ro',ecolor='b', lw=2, capsize=2, capthick=1)
    plt.xticks(x_axis,algo_label,fontsize=18)
    ax.tick_params(axis="y", labelsize=18)
    # ax.legend(loc='upper right')
    #add x and y labels
    ax.set_ylim(0,10)

    ax.set_ylabel('Relative Error Percentage (%)',fontsize=20)
    ax.set_xlabel('Real world graph datasets',fontsize=20)

    fig.suptitle('Accuracy of TETRIS',fontsize=20)

    #show plot
    # plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # fig.set_size_inches(6,6)
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
        out_filename = "-accuracy-"
        data_dir.append("output/plot_data/accuracy_plot_data/EdgeEstimated/"+file+".edges/")
    plot_accuracy(file_name, data_dir, out_filename)