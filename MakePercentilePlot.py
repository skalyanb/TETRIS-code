import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt

dirname = os.path.dirname
sys.path.append(dirname(dirname(os.path.realpath(__file__))))

def skip_rows (index):
    if index <4 or index > 5:
        return True
    return False

def plot_estimates (data_dir, out_filename, x_max):
    X=[]
    Y=[]
    max_err_percent = 10
    min_err_percent = 10

    max_band_percent = 1.5
    min_band_percent = 1.5

    for i, file in enumerate(os.listdir(data_dir)):
        path = os.path.join(data_dir, file)
        if os.path.isdir(path):
            # skip directories
            continue
        if not file.startswith('.'):
            data_file = data_dir + file
            df_tri = pd.read_csv(data_file, sep=",",skiprows=lambda x: skip_rows(x))
            triangle_count = df_tri["triangles"][0]
            df_rawdata = pd.read_csv(data_file, sep=",",skiprows=15)
            # keep only the data that are withing 80 percentile
            upper = df_rawdata["triangle_estimate"].quantile(.75)
            lower = df_rawdata["triangle_estimate"].quantile(.25)
            df_rawdata = df_rawdata[df_rawdata["triangle_estimate"] < upper]
            df_rawdata = df_rawdata[df_rawdata["triangle_estimate"] > lower]
            Y.extend(df_rawdata["triangle_estimate"])
            X.extend(df_rawdata["fraction_of_edges_seen"])

    fig, ax = plt.subplots()

    y_band_min = triangle_count*(1-min_band_percent/100.0)
    y_band_max = triangle_count*(1+max_band_percent/100.0)
    ax.axhline(y=triangle_count, linewidth=2, color='r', linestyle='-',label='Count')
    ax.axhline(y=y_band_min, color='r', linestyle='--')
    ax.axhline(y=y_band_max, color='r', linestyle='--')

    #scatter plot
    ax.scatter(X, Y, s=20, c='blue', marker='o')

    #change axes ranges
    y_min = triangle_count * (1-min_err_percent/100.0)
    y_max = triangle_count * (1+max_err_percent/100.0)
    ax.set_ylim(y_min, y_max)

    x_min = 0.0
    x_max = x_max
    ax.set_xlim(x_min,x_max)

    relative_pos_band_exact = 1.0*min_err_percent / (min_err_percent + max_err_percent)
    relative_pos_band_min = 1.0* (min_err_percent-min_band_percent) / (min_err_percent + max_err_percent)
    relative_pos_band_max = 1.0* (min_err_percent+max_band_percent) / (min_err_percent + max_err_percent)

    ax.text(1.02, relative_pos_band_exact, '0%',
            horizontalalignment='left',
            verticalalignment='center',
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, relative_pos_band_max, '+'+str(max_band_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, relative_pos_band_min, '-'+str(min_band_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, 1, '+'+str(max_err_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, 0, '-'+str(min_err_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            #rotation='vertical',
            transform=ax.transAxes
            )


    #add title
    ax.set_title('Accuracy vs Observed Graph')

    #add x and y labels
    ax.set_xlabel('Percentage of Edges Visited')
    ax.set_ylabel('Triangle Estimates')

    #show plot
    #plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")

    fig.savefig("output/plots/variance/small_fraction/"+out_filename+"-"+timestr+".eps",format='eps')




if __name__ == "__main__":

    # file_name = "soc-flickr-und"
    # x_max = 1.1


    file_name = "soc-sinaweibo"
    x_max = 1.1

    # file_name = "soc-flickr"
    # x_max = 1.1

    # file_name = "soc-twitter-konect"
    # x_max = 0.11

    # file_name = "soc-livejournal"
    # x_max = 1.1

    # file_name = "soc-orkut"
    # x_max = 1.1

    # file_name = "socfb-A-anon"
    # x_max = 1.1

    data_dir = "output/plot_data/variance_plot_data/"+file_name+".edges/EstTriByRWandWghtedSampling/small_fraction/"
    out_filename = file_name

    plot_estimates(data_dir,out_filename,x_max)

    #plot_comparison(data_dir_1,data_dir_2)

