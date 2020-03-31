#! bin/python3
import os
import csv
import pandas as pd 
import seaborn as sns 
import numpy as np
import matplotlib.pyplot as plt
import gc
import pickle
import time
from scipy.stats import pearsonr, spearmanr
from utils import *
gc.disable()

def main():
    not_found = []
    cluster_vals = []
    spearmanr_vals = []
    pearsonr_vals = []
    abc_scores_vals ={}
    distance_vals = {}
    dist_not_found_1 = []
    clusters = pd.read_csv("non_expressed/sorted_nonexpressed_clusters.txt", sep="\t", header=None)
    cells = pd.read_csv("SeventyOne_CellTypes.txt", sep="\t", header=None)
    save_outdir = "/srv/scratch/kmualim/non_expressed_links/"
    for name in clusters[0]:
        print("Cluster", name)
        outdir="non_expressed"
        try:
            dat = pd.read_csv("non_expressed/{}_Intersected.tsv.gz".format(name), sep="\t", compression="gzip")
            data = dat.rename(columns={'DHS.RPM.TSS1Kb':'Gene.DHS.RPM.TSS1Kb', 'H3K27ac.RPM.TSS1Kb':'Gene.H3K27ac.RPM.TSS1Kb'})
            target_gene = data['TargetGene'].drop_duplicates().values[0]
            pearsonr,spearmanr, abc_scores, distance = makeDataFrame(data, cells[0], name)
            cluster_vals.append(name)
            pearsonr_vals.append(pearsonr)
            spearmanr_vals.append(spearmanr)
            abc_scores_vals[name] = abc_scores
            distance_vals[name] = distance
        except:
            not_found.append(name)
            print("Could not find Cluster {}".format(name))
            continue
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
    
    with open(os.path.join(save_outdir, "clusters_notfound_cellines_{}.txt".format(name)), "w") as f:
        for listitem in not_found:
            f.write('{}\n'.format(str(listitem)))
        f.close()

    print("Saving ABCScoreValues")
    pickle.dump(abc_scores_vals, open(os.path.join(save_outdir, "ABCScoreValues_Cluster_{}.p".format(name)), "wb"), protocol=-1)
    print("Saving DistanceValues")
    pickle.dump(distance_vals, open(os.path.join(save_outdir, "DistanceValues_Cluster_{}.p".format(name)), "wb"), protocol=-1)

    print(len(spearmanr_vals))
    print(len(pearsonr_vals))
    print(len(cluster_vals))
    print(len(abc_scores_vals))
    print(len(distance_vals))
    
    # plot ABC behaviour across K562
    plot_outdir = "/oak/stanford/groups/akundaje/projects/ABC_links/plots/CellTypes_EnhancerPromoterPropertiesPlots/Version2_plots/non_expressed_links/plots"
    celltypes = pd.read_csv("SixtyEight_CellTypes.txt", sep="\t", header=None)
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
    
def grabExpressionValues(cells, target_gene, expression_files):
    values = []
#    import pdb; pdb.set_trace()
    for i,j in zip(cells,expression_files):
        try:
            express_val = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/scripts/cluster_expression_files/{}.{}.TPM.tsv".format(i, j), sep="\t", header=None)
        except:
            express_val = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/scripts/cluster_expression_files/{}.{}.TPM.tsv".format(i.split("_")[1], j), sep="\t", header=None)
        val = express_val.loc[express_val[0]==target_gene]
        if len(val)!=0:
            val_ = np.log(val.values[0][1]+0.001)
        else:
            val_ = np.log(0.0+0.001)
        values.append(val_)
    return values 

def grabValuesFromAppropriateColumn(data, column, cells):
    values = []
    # grab values for K562 first 
    k562 = data[data[column].notnull()]
    val = np.log(k562[column].values[0]+0.001)
    values.append(val)

    for i in cells: 
        x = data[data[column+'_{}'.format(i)].notnull()]
        if len(x) != 0:
            val = np.log(x[column+'_{}'.format(i)].values[0]+0.001)
        else:
            val = 0
        values.append(val)
    return values

def calculateMetrics(list_metrics):
    import itertools
    pearson_vals = []
    spearman_vals = []
    
    for a, b in itertools.combinations(list_metrics, 2):
        r, _ = pearsonr(a, b)
        x, _ = spearmanr(a,b)
        pearson_vals.append(r)
        spearman_vals.append(x)
    return pearson_vals, spearman_vals

def corrfunc(x,y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = pearsonr(x, y)
    x, _ = spearmanr(x, y)
    ax = ax or plt.gca()
    # Unicode for lowercase rho (ρ)
    rho = '\u03C1'
    ax.annotate(f'p {rho} = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)
    ax.annotate(f's {rho} = {x:.2f}', xy=(.1, 1.0), xycoords=ax.transAxes)

def get_spearmanr(x1, y):
    r, _ = pearsonr(x1, y)
    x, _ = spearmanr(x1, y)
    return r,x

def getActivity(data, cells):
    # activity for K562
    data['Promoter_activity_base'] = np.sqrt(data['Gene.DHS.RPM.TSS1Kb']*data['Gene.H3K27ac.RPM.TSS1Kb'])
    # activity for cells 
    for i in cells:
        data['Promoter_activity_base_{}'.format(i)] = np.sqrt(data['Gene.DHS.RPM.TSS1Kb_{}'.format(i)]*data['Gene.H3K27ac.RPM.TSS1Kb_{}'.format(i)])
    return data

def getBinCategory(distance):
    categories = [] 
    for i in distance:
        if i<10000 and i>0:
            categories.append('0KbTo10Kb')
        elif i>10000 and i<100000:
            categories.append('10KbTo100Kb')
        elif i>100000 and i<1000000:
            categories.append('100KbTo1000Kb')
        else:
            categories.append('GreaterThan1000Kb')
    assert len(categories) == len(distance)
    return categories 

def makeDataFrame(data, cells, name):
    dist_not_found=[]
    print("Making dataframe")
    df = pd.DataFrame()
    dhs= grabValuesFromAppropriateColumn(data,'normalized_dhs', cells)
    h3k27ac= grabValuesFromAppropriateColumn(data,'normalized_h3K27ac', cells)
    dhs_tss= grabValuesFromAppropriateColumn(data,'Gene.DHS.RPM.TSS1Kb', cells)
    h3k27ac_tss= grabValuesFromAppropriateColumn(data,'Gene.H3K27ac.RPM.TSS1Kb', cells)    
    target_gene = data['TargetGene'].drop_duplicates().values[0] 
    data = getActivity(data, cells)
    promoter_activity= grabValuesFromAppropriateColumn(data,'Promoter_activity_base', cells)
    activity= grabValuesFromAppropriateColumn(data,'activity_base', cells)
    abc_scores = grabValuesFromAppropriateColumn(data, 'ABC.Score', cells)
    distance = grabValuesFromAppropriateColumn(data,'distance', cells)
    pearsonr_vals, spearmanr_vals = get_spearmanr(promoter_activity, activity)

    files_outdir="/oak/stanford/groups/akundaje/projects/ABC_links/plots/CellTypes_EnhancerPromoterPropertiesPlots/Version2_plots/non_expressed_links/plots"
    outdir = "/srv/scratch/kmualim/non_expressed_links/"
    data= {'logEnhancers.DHS.RPM': dhs,
            'logEnhancers.H3K27ac.RPM':h3k27ac,
            'logDHS.RPM.TSS1kb':dhs_tss,
            'logH3K27ac.RPM.TSS1kb':h3k27ac_tss,
            'activity_base': activity,
            'promoter_activity': promoter_activity,
            'logABC.Score': abc_scores}
            #'Pearson Correlation': pearsonr_vals,
            #'Spearman Correlation': spearmanr_vals}
    df = pd.DataFrame(data)
    df['CellTypes'] = ['K562']+['other']*len(cells)
    print("Plotting the Pairplot")
    list_metrics = [dhs, h3k27ac, dhs_tss, h3k27ac_tss, promoter_activity]
    pearsonr, spearmanr= calculateMetrics(list_metrics)
    g = sns.pairplot(apply_jitter(df), x_vars={'logEnhancers.DHS.RPM', 'logEnhancers.H3K27ac.RPM', 'activity_base'}, y_vars={'logDHS.RPM.TSS1kb', 'logH3K27ac.RPM.TSS1kb', 'promoter_activity'}, hue='CellTypes', diag_kind='hist')
    g.savefig(os.path.join(files_outdir, "{}_{}.Correlation.png".format(name, target_gene)))
    g.savefig(os.path.join(outdir, "{}_{}.Correlation.pdf".format(name, target_gene)), format='pdf')
    plt.close('all')
    with open(os.path.join(outdir, "{}_{}_celllines_correlation_values.txt".format(name, target_gene)), "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(pearsonr, spearmanr))
    print("Saving Plots...") 
    df.to_csv(os.path.join(outdir, "{}_{}_celllines.tsv.gz".format(name, target_gene)), compression="gzip", sep="\t", index=False)
    
    return pearsonr_vals, spearmanr_vals, abc_scores, distance#, dist_not_found
if __name__=='__main__':
   main() 

