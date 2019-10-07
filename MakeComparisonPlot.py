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



def plot_comparison(data_dir, filename, title_info):

    X = []
    Y= []

    for i,dir in enumerate(data_dir):
        tempX, tempY = populate_data(data_dir)
        X.append(tempX)
        Y.append(tempY)
    # X_1, Y_1 = populate_data(data_dir_1)
    #
    # X_2, Y_2 = populate_data(data_dir_2)
    #
    # X_3, Y_3 = populate_data(data_dir_3)
    #
    # X_4, Y_4 = populate_data(data_dir_4)

    #scatter plot
    fig, ax = plt.subplots()
    ax.plot(X[0], Y[0], c='blue', marker='s', label='Our Algo')

    ax.plot(X[1], Y[1], c='red', marker='o', label='UES-&-Count')

    ax.plot(X[2], Y[2], c='red', marker='^', label='RW-&-Count')

    ax.plot(X[3], Y[3], c='green', marker='D', label='UES-sparsify')

    ax.set_xlim(2,20)
    ax.legend(loc='upper right')

    # plt.ylim(y_min, y_max)

    #add title
    ax.set_title('Comparison against Baselines\n'+title_info)

    #add x and y labels
    ax.set_xlabel('Percentage of Edges Visited')
    ax.set_ylabel('Median Error % (100 runs)')

    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    #fig.savefig("output/plots/comparison/"+filename+timestr+".eps",format='eps')

if __name__ == "__main__":


    file_name = []
    title_info = []

    # f_name = "soc-flickr-und"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 16M edges, 1.7M vertices")

    # f_name = "socfb-A-anon"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 24M edges, 3M vertices")

    f_name = "soc-orkut"
    file_name.append(f_name)
    title_info.append(f_name + ": 106M edges, 3M vertices")

    # f_name = "soc-sinaweibo"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 260M edges, 58M vertices")
    #
    # f_name = "soc-twitter-konect"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.2B edges, 41M vertices")
    #
    # f_name = "soc-friendster"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")

    algo_list = []

    algo_list.append("EstTriByRWandWghtedSampling")
    algo_list.append("EstTriByEdgeSampleAndCount")
    algo_list.append("EstTriByRWAndCountPerEdge")
    algo_list.append("EstTriByUniformSampling")

    for i,file in enumerate(file_name):
        data_dir = []
        out_filename = file + "-comparison-"
        for j,algo in enumerate(algo_list):
            data_dir.appned("output/plot_data/variance_plot_data/"+file+".edges/"+algo+"/fixed_seed/")
        plot_comparison(data_dir, out_filename, title_info[i])