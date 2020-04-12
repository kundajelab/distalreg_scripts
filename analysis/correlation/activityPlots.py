#! bin/python3
import os, sys
import csv
import tqdm
import pandas as pd 
import seaborn as sns 
import numpy as np
import matplotlib.pyplot as plt
import gc
from multiprocessing import Pool
import pickle
import time
from scale_values import scale_values
import argparse
from tools import *
from scipy.stats import pearsonr, spearmanr

pd.options.mode.chained_assignment = None
gc.disable()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--clusters', help="File with enhancer regions")
    parser.add_argument('--save_outdir', help="Path to save files")
    parser.add_argument('--celltypes', help="File with list of celltypes")
    args = parser.parse_args()
    return args 


def main():
    args = parse_args()
#    clusters = "../apr2_datafiles/expressed_files_apr2/expressed_clusters.txt"
#    clusters_1 = pd.read_csv(clusters, sep="\t", header=None)
    clusters = sys.argv[1]
    clusters_1 = pd.read_csv(clusters, sep="\t", header=None)
    with Pool(10) as p:
        p.map(generateIndividualClusterFiles, list(clusters_1[0]))

def generateIndividualClusterFiles(clusters):
    files_outdir="/oak/stanford/groups/akundaje/projects/ABC_links/plots/CellTypes_EnhancerPromoterPropertiesPlots/Version2_plots/expressed_links/plots"
    plot_outdir = "/oak/stanford/groups/akundaje/projects/ABC_links/plots/CellTypes_EnhancerPromoterPropertiesPlots/Version2_plots/expressed_links/plots"
    not_found = []
    cluster_vals = []
    spearmanr_vals = []
    pearsonr_vals = []
    abc_scores_vals ={}
    distance_vals = {}
    dist_not_found_1 = []
#    clusters = pd.read_csv(args.clusters, sep="\t", header=None)
    celltypes = pd.read_csv("cells.txt", sep="\t", header=None)
    save_outdir = "/srv/scratch/kmualim/expressed_links/accumulate"
#    for name in clusters:
    name = clusters
    print("Cluster", name)
    datafile = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/apr2_datafiles/expressed_files_apr2/expressed_cluster_files_apr3/{}".format(name), sep="\t",compression="gzip")
    # filtering dataframe
    index = datafile[datafile['TargetGene'] == 'TargetGene'].index
    datafile.drop(index, inplace=True)
    datafile['activity_base'] = datafile['activity_base'].astype('float')
    dat = datafile.loc[datafile['activity_base']>0.0]
    if len(dat) > 0:
        data = dat.rename(columns={'Gene.DHS.RPM.TSS1Kb':'Gene.DHS.RPM.TSS1Kb_osteoblast', 'Gene.H3K27ac.RPM.TSS1Kb':'Gene.H3K27ac.RPM.TSS1Kb_osteoblast', 'DHS.RPM.TSS1Kb':'Gene.DHS.RPM.TSS1Kb', 'H3K27ac.RPM.TSS1Kb':'Gene.H3K27ac.RPM.TSS1Kb'})
        # filter rows with no celltype
        data_1 = data[data['cellType'].notna()]
        # filter for specific target genes
        unique_data = data_1['TargetGene'].drop_duplicates() 
        for targetgene in unique_data:
            print(targetgene)
            data_2 = data_1.drop_duplicates()
            subset_data = data_2.loc[data_2['TargetGene']==targetgene]
            pearsonr,spearmanr, abc_scores, distance = makeDataFrame(subset_data, celltypes[0], name, files_outdir)
            cluster_vals.append(name)
            pearsonr_vals.append(pearsonr)
            spearmanr_vals.append(spearmanr)
            abc_scores_vals[name] = abc_scores
            distance_vals[name] = distance
    else:
        not_found.append(name)
        print("Could not find/Not Expressed Cluster {}".format(name))
            
    print("Saving Summary")
    with open(os.path.join(save_outdir, "Summary_{}.txt".format(name)), "w") as f:
        for spearmanr, pearsonr, cluster in zip(spearmanr_vals, pearsonr_vals, cluster_vals):
            f.write(str(cluster))
            f.write("\t")
            f.write(str(pearsonr))
            f.write("\t")
            f.write(str(spearmanr))
            f.write("\n")
        f.close()

    print("Clusters that were not found")
    with open(os.path.join(save_outdir, "clusters_notfound_cellines_{}.txt".format(name)), "w") as f:
        for listitem in not_found:
            f.write('{}\n'.format(str(listitem)))
        f.close()

    print("Saving ABCScoreValues")
    pickle.dump(abc_scores_vals, open(os.path.join(save_outdir, "ABCScoreValues_Cluster_{}.p".format(name)), "wb"), protocol=-1)
    print("Saving DistanceValues")
    pickle.dump(distance_vals, open(os.path.join(save_outdir, "DistanceValues_Cluster_{}.p".format(name)), "wb"), protocol=-1)

    
    for cell in range(len(celltypes)):
        k562_abc_scores = [score[cell] for score in abc_scores_vals.values()]
        k562_dist_scores = [distance[cell] for  distance in distance_vals.values()]
        dist_cat = getBinCategory(k562_dist_scores)

        # abc scores
        df = pd.DataFrame()
        df['logABC.score'] = np.concatenate((k562_abc_scores, k562_abc_scores))
        df['Correlation'] = np.concatenate((spearmanr_vals, pearsonr_vals))
        df['Value'] = ['Spearman']*len(spearmanr_vals) + ['Pearson']*len(pearsonr_vals)
        g = sns.jointplot(x='logABC.score', y='Correlation', col='Value', data=df)
        g.savefig(os.path.join(plot_outdir, "logABC.Score.Correlation.{}.png".format(celltypes[cell])))
        g.savefig(os.path.join(save_outdir, "logABC.Score.Correlation.{}.pdf".format(celltypes[cell])), format='pdf')
        plt.close('all')
    
        df = pd.DataFrame()
        df['DistanceOfEnhancerFromPromoter'] = np.concatenate((dist_cat, dist_cat))
        df['Correlation'] = np.concatenate((spearmanr_vals, pearsonr_vals))
        df['Value'] = ['Spearman']*len(spearmanr_vals) + ['Pearson']*len(pearsonr_vals)
        g = sns.catplot(x='DistanceOfEnhancerFromPromoter', y='Correlation', col='Value', data=df, kind='violin')
        g.savefig(os.path.join(plot_outdir, "DistanceOfEnhancerFromPromoter.Correlation.{}.png".format(celltypes[cell])))
        g.savefig(os.path.join(save_outdir, "DistanceOfEnhancerFromPromoter.Correlation.{}.pdf".format(celltypes[cell])), format='pdf')
        plt.close('all')
    

if __name__=='__main__':
   main() 

