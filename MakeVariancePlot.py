import os
import sys
import time
import pandas as pd
import matplotlib.pyplot as plt
import statistics


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

def skip_rows (index):
    if index <4 or index > 5:
        return True
    return False

def plot_estimates (data_dir, out_filename, x_max, title_info):
    edges_seen=[]
    tri_est=[]
    max_err_percent = 20
    min_err_percent = 20

    max_band_percent = 5
    min_band_percent = 5

    triangle_count = 0

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
            tri_est.extend(df_rawdata["triangle_estimate"])
            edges_seen.extend(df_rawdata["fraction_of_edges_seen"])

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom = 0.19, top = 0.8, right = 0.85, left = 0.15, hspace = 0.38)

    y_band_min = triangle_count*(1-min_band_percent/100.0)
    y_band_max = triangle_count*(1+max_band_percent/100.0)
    ax.axhline(y=triangle_count, linewidth=2, color='r', linestyle='-',label='Count')
    ax.axhline(y=y_band_min, color='r', linestyle='--')
    ax.axhline(y=y_band_max, color='r', linestyle='--')

    #scatter plot
    ax.scatter(edges_seen, tri_est, s=20, c='blue', marker='o')

    #change axes ranges
    y_min = triangle_count * (1-min_err_percent/100.0)
    y_max = triangle_count * (1+max_err_percent/100.0)
    ax.set_ylim(y_min, y_max)

    x_min = 0.1
    x_max = 3.2
    ax.set_xlim(x_min,x_max)

    relative_pos_band_exact = 1.0*min_err_percent / (min_err_percent + max_err_percent)
    relative_pos_band_min = 1.0* (min_err_percent-min_band_percent) / (min_err_percent + max_err_percent)
    relative_pos_band_max = 1.0* (min_err_percent+max_band_percent) / (min_err_percent + max_err_percent)

    ax.text(1.02, relative_pos_band_exact, '0%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, relative_pos_band_max, '+'+str(max_band_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, relative_pos_band_min, '-'+str(min_band_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, 1, '+'+str(max_err_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )
    ax.text(1.02, 0, '-'+str(min_err_percent)+'%',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=20,
            #rotation='vertical',
            transform=ax.transAxes
            )


    #add title
    ax.set_title(title_info,fontsize=20)

    #add x and y labels
    ax.set_xlabel('Percentage of queries %',fontsize=20)
    ax.set_ylabel('Triangle Estimate',fontsize=20)

    ax.tick_params(axis="y", labelsize=20)
    ax.tick_params(axis="x", labelsize=20)
    #show plot
    plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")

    fig.savefig("output/plots/variance/"+out_filename+"-"+timestr+".eps",bbox_inches="tight")




if __name__ == "__main__":

    x_max = 3.1

    file_name = []
    title_info = []

    # f_name = "soc-flickr-und"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 32M edges")

    f_name = "soc-orkut"
    file_name.append(f_name)
    title_info.append(f_name + ": 213M edges")

    f_name = "soc-sinaweibo"
    file_name.append(f_name)
    title_info.append(f_name + ": 523M edges")
    #
    # f_name = "soc-twitter-konect"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 2.4B edges")

    # f_name = "socfb-A-anon"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 24M edges, 3M vertices")
    # #

    # f_name = "soc-friendster"
    # file_name.append(f_name)
    # title_info.append(f_name + ": 1.8B edges, 65M vertices")
    # #
    for i,file in enumerate(file_name):
        data_dir = "output/plot_data/baseline_TETRIS_SRW_MCMC/"+file+".edges.csr/TETRIS/"
        out_filename = file

        plot_estimates(data_dir,out_filename,x_max,title_info[i])

    #plot_comparison(data_dir_1,data_dir_2)

