#! bin/python3 

import pandas as pd 
import numpy as np
import math
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

clusters = pd.read_csv("sorted_cluster_1.txt", sep="\t", header=None)
genes_seen = []

def scale_values(single_data):

    # iterating through columns that contain values
    for col in single_data.columns[:7]:
        max_val = np.max(single_data[col])
        min_val = np.min(single_data[col])
        single_data[col+"_normalized"] = [(row-min_val)/(max_val-min_val) for row in single_data[col]]
        median = single_data[col+"_normalized"].median()
        single_data[col+"_median_centered"] = single_data[col+"_normalized"] - median
    return single_data

def process(single_data, outdir, file_, name):
    single_data = scale_values(single_data)
    single_data.to_csv(os.path.join(outdir, "{}_ScaledForNull.tsv.gz".format(file_)), compression="gzip", sep="\t", index=False)
    ax = sns.pairplot(x_vars={'logEnhancers.DHS.RPM_median_centered', 'logEnhancers.H3K27ac.RPM_median_centered', 'activity_base_median_centered'}, y_vars={'logDHS.RPM.TSS1kb_median_centered', 'logH3K27ac.RPM.TSS1kb_median_centered', 'promoter_activity_median_centered'}, data=single_data, hue='CellTypes')
    def plot_unity(xdata, ydata, **kwargs):
        mn = min(xdata.min(), ydata.min())
        mx = max(xdata.max(), ydata.max())
        points = np.linspace(mn, mx, 100)
        plt.gca().plot(points, points, color='k', marker=None,
                    linestyle='--', linewidth=1.0)
    ax.map_offdiag(plot_unity)
    ax.savefig(os.path.join(outdir,"{}_ScaledAndMedianCenteredPlot.png".format(name)))
    plt.close("all")
        #return single_data


def main():
    no_clusters=[]
    for cluster in clusters[0]:
        print("Running Gene: ",cluster)
        outdir = "/oak/stanford/groups/akundaje/projects/ABC_links/plots/CellTypes_EnhancerPromoterPropertiesPlots/Version2_plots/thresholded_plots/"
        try:
            files = glob.glob(os.path.join(outdir, "{}*.tsv.gz".format(cluster)))
            single_data = pd.read_csv(files[0], compression="gzip", sep="\t")
            genes_seen.append(cluster)
            process(single_data, outdir, files[0], str(files[0]).split(".")[0])
        except:
            no_clusters.append(cluster)
        
        with open("cluster_not_found.txt", "w") as f:
            for i in no_clusters:
                f.write(str(i))
                f.write("\n")
            f.close()

if __name__=="__main__":    
    main(clusters)
