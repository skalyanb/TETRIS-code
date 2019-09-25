import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
dirname = os.path.dirname
sys.path.append(dirname(dirname(os.path.realpath(__file__))))

def skip_rows (index):
    if index <4 or index > 5:
        return True
    return False

def plot_estimates (data_dir):
    X=[]
    Y=[]
    for i, file in enumerate(os.listdir(data_dir)):
        data_file = data_dir + file
        df_tri = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x))
        triangle_count = df_tri["triangles"][0]
        df_rawdata = pd.read_csv(data_file, sep=",",skiprows=15)
        Y.extend(df_rawdata["triangle_estimate"])
        X.extend(df_rawdata["fraction_of_edges_seen"])


    plt.axhline(y=triangle_count, color='r', linestyle='-')
    plt.axhline(y=triangle_count*0.9, color='r', linestyle='-')
    plt.axhline(y=triangle_count*1.1, color='r', linestyle='-')

    #scatter plot
    plt.scatter(X, Y, s=20, c='blue', marker='o')

    #change axes ranges
    y_min = triangle_count * 0.7
    y_max = triangle_count * 1.3
    plt.xlim(0,30)
    plt.ylim(y_min, y_max)

    #add title
    plt.title('Accuracy vs Visits')

    #add x and y labels
    plt.xlabel('Percentage of Edges Visited')
    plt.ylabel('Triangle Estimates')

    #show plot
    plt.show()


if __name__ == "__main__":
    DEFAULT = "output/soc-flickr.edges"

    data_dir = DEFAULT if len(sys.argv) == 1 else sys.argv[1]
    data_dir += "/"

    #fig, ax = plt.subplots(1, 1)

    # key, ylabel = plot_estimates(ax, data_dir)
    plot_estimates(data_dir)

    # fig.legend(handles=key)
    # plt.ylabel(ylabel)
    # plt.show()