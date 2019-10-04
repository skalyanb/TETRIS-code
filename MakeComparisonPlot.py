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

    df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
    edges_seen = df_edges_seen.iloc[0,0]

    return median_error,edges_seen

def populate_data (data_dir):
    X = []
    Y=[]
    for i, file in enumerate(os.listdir(data_dir)):
        data_file = data_dir + file
        if valid_file(data_dir,file):
            median_error, edges_seen = parse_outfile(data_file)
            Y.append(median_error)
            X.append(edges_seen)

    return X,Y



def plot_comparison(data_dir_1,data_dir_2,data_dir_3,data_dir_4, filename, title_info):

    X_1, Y_1 = populate_data(data_dir_1)

    X_2, Y_2 = populate_data(data_dir_2)

    X_3, Y_3 = populate_data(data_dir_3)

    X_4, Y_4 = populate_data(data_dir_4)

    #scatter plot
    fig, ax = plt.subplots()
    ax.plot(X_1, Y_1, c='blue', marker='s', label='Our Algo')

    ax.plot(X_2, Y_2, c='red', marker='o', label='UES-&-Count')

    ax.plot(X_3, Y_3, c='red', marker='^', label='RW-&-Count')

    ax.plot(X_4, Y_4, c='green', marker='D', label='UES-sparsify')

    ax.set_xlim(2,20)
    ax.legend(loc='upper right')

    # plt.ylim(y_min, y_max)

    #add title
    ax.set_title('Comparison against Baselines\n'+title_info)

    #add x and y labels
    ax.set_xlabel('Percentage of Edges Visited')
    ax.set_ylabel('Median Error % (100 runs)')

    #show plot
    # plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps')

if __name__ == "__main__":

    # file_name = "soc-flickr"
    # title_info = "soc-flickr: 6.4M edges, 514K vertices"
    # x_min = 2
    # x_max = 20

    file_name = "soc-flickr-und"
    title_info = "soc-flickr: 31M edges, 1.7M vertices"
    x_min = 2
    x_max = 20


    data_dir_1 = "output/plot_data/comparison_plot_data/"+file_name +".edges/EstTriByRWandWghtedSampling/"
    data_dir_2 = "output/plot_data/comparison_plot_data/"+file_name+".edges/EstTriByEdgeSampleAndCount/"
    data_dir_3 = "output/plot_data/comparison_plot_data/"+file_name +".edges/EstTriByRWAndCountPerEdge/"
    data_dir_4 = "output/plot_data/comparison_plot_data/"+file_name +".edges/EstTriByUniformSampling/"
    filename = file_name + "-comparison-"


    plot_comparison(data_dir_1,data_dir_2,data_dir_3, data_dir_4, filename, title_info)

