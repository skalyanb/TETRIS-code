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

def skip_rows_err (index):
    if index <10 or index > 10:
        return True
    return False


def skip_rows_seen (index):
    if index <13 or index > 13:
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
    fig = plt.figure()
#    plt.show()
    fig.save(data_dir+".png")

def plot_comparison(data_dir_1,data_dir_2):
    X_1=[]
    Y_1=[]
    X_2=[]
    Y_2=[]

    for i, file in enumerate(os.listdir(data_dir_1)):
        data_file = data_dir_1 + file
        print (data_file)

    #df_tri_count = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x),header=None, usecols=[2])
        #triangle_count = df_tri_count[0][0]

        df_median_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1])
        median_error = df_median_error.iloc[0,0]

        df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
        edges_seen = df_edges_seen.iloc[0,0]

        Y_1.append(median_error)
        X_1.append(edges_seen)

    for i, file in enumerate(os.listdir(data_dir_2)):
        data_file = data_dir_2 + file
        print (data_file)

        #df_tri_count = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x),header=None, usecols=[2])
        #triangle_count = df_tri_count[0][0]

        df_median_error = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_err(x),header=None, usecols=[1])
        median_error = df_median_error.iloc[0,0]

        df_edges_seen = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows_seen(x),header=None, usecols=[0])
        edges_seen = df_edges_seen.iloc[0,0]

        Y_2.append(median_error)
        X_2.append(edges_seen)

        #scatter plot
    plt.scatter(X_1, Y_1, s=20, c='blue', marker='o')

    plt.scatter(X_2, Y_2, s=20, c='red', marker='o')

    #change axes ranges
    # y_min = triangle_count * 0.7
    # y_max = triangle_count * 1.3
    plt.xlim(0,40)
    # plt.ylim(y_min, y_max)

    #add title
    plt.title('Accuracy vs Visits')

    #add x and y labels
    plt.xlabel('Percentage of Edges Visited')
    plt.ylabel('Triangle Estimates')

    #show plot
    #fig = plt.figure()
    plt.show()
    # fig.save(data_dir+".png")


if __name__ == "__main__":
    data_dir_1 = "output/socfb-A-anon.edges/EstTriByRWWfgtdSamp/"
    #data_dir_2 = "output/socfb-A-anon.edges/EstTriByRWAndCount/"
    data_dir_2 = "output/socfb-A-anon.edges/EstTriByRWAndNeighborSample/"

    #data_dir = DEFAULT if len(sys.argv) == 1 else sys.argv[1]
    # if len(sys.argv) == 2:
    #     data_dir_2 = sys.argv[2]
    #     data_dir_2 += "/"
    # data_dir += "/"

    #fig, ax = plt.subplots(1, 1)

    # key, ylabel = plot_estimates(ax, data_dir)
    #plot_estimates(data_dir)

    plot_comparison(data_dir_1,data_dir_2)

    # fig.legend(handles=key)
    # plt.ylabel(ylabel)
    # plt.show()