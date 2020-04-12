import numpy as np
import os
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
from scipy.stats import pearsonr, spearmanr


def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def apply_jitter(df):
    for i in range(len(df.columns)-1):
        df.iloc[:, i] = rand_jitter(df.iloc[:, i])
    return df

    
def grabExpressionValues(cells, target_gene, expression_files):
    values = []
    for i,j in zip(cells,expression_files):
        try:
            express_val = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/scripts/cluster_expression_files/{}.{}.TPM.tsv".format(i, j), sep="\t", header=None)
        except:
            express_val = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/scripts/cluster_expression_files/{}.{}.TPM.tsv".format(i.split("_")[1], j), sep="\t", header=None)
        val = express_val.loc[express_val[0]==target_gene]
        if len(val)!=0:
            val_ = np.arcsinh(val.values[0][1]+0.001)
        else:
            val_ = np.arcsinh(0.0+0.001)
        values.append(val_)

def grabValuesFromColumns(data, column, cells):
    values = []
    total_col = data.columns[data.columns.str.contains(column+"_")][0]
    k562_data = data.fillna(0.0)
    val = np.arcsinh(np.float(k562_data[column].values[0])+0.001)
    values.append(val)
     
    for i in cells:
        matched = data.loc[data['cellType']==i]
        if len(matched) != 0:
            val = np.arcsinh(np.float(matched[total_col].values[0])+0.001)
        else:
            val = np.arcsinh(0.0)
        values.append(val)
    return values
        

def grabValuesFromAppropriateColumn(data, column, cells):
    values = []
    # grab values for K562 first 
    k562 = data[data[column].notnull()]
    val = np.log(np.float(k562[column].values[0])+0.001)
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

def getActivityColumn(data, cells):
    data['Promoter_activity_base'] = np.sqrt(data['Gene.DHS.RPM.TSS1Kb'].astype('float')*data['Gene.H3K27ac.RPM.TSS1Kb'].astype('float'))
    column = data.columns[data.columns.str.contains('Gene.DHS.RPM.TSS1Kb_')][0]
    h3k27ac_col = data.columns[data.columns.str.contains('Gene.H3K27ac.RPM.TSS1Kb_')][0]
    for i in cells:
        matched = data.loc[data['cellType']==i]
        if len(matched) != 0:
            data['Promoter_activity_base_{}'.format(i)] = np.sqrt(np.float(matched[column].values[0])*np.float(matched[h3k27ac_col].values[0]))
        else:
            data['Promoter_activity_base_{}'.format(i)] = 0.0
    return data 

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

def makeDataFrame(data, cells, name, files_outdir):
    dist_not_found=[]
    print("Making dataframe")
    df = pd.DataFrame()
    dhs= grabValuesFromColumns(data,'normalized_dhs', cells)
    h3k27ac= grabValuesFromColumns(data,'normalized_h3K27ac', cells)
    dhs_tss= grabValuesFromColumns(data,'Gene.DHS.RPM.TSS1Kb', cells)
    h3k27ac_tss= grabValuesFromColumns(data,'Gene.H3K27ac.RPM.TSS1Kb', cells)
    target_gene = data['TargetGene'].drop_duplicates().values[0] 
    data = getActivityColumn(data, cells)
    promoter_activity= grabValuesFromColumns(data,'Promoter_activity_base', cells)
    activity= grabValuesFromColumns(data,'activity_base', cells)
    abc_scores = grabValuesFromColumns(data, 'ABC.Score', cells)
    distance = grabValuesFromColumns(data,'distance', cells)
    pearsonr_vals, spearmanr_vals = get_spearmanr(promoter_activity, activity)

    outdir = "/srv/scratch/kmualim/expressed_links/accumulate"
    data= {'logEnhancers.DHS.RPM': dhs,
            'logEnhancers.H3K27ac.RPM':h3k27ac,
            'logDHS.RPM.TSS1kb':dhs_tss,
            'logH3K27ac.RPM.TSS1kb':h3k27ac_tss,
            'activity_base': activity,
            'promoter_activity': promoter_activity,
            'logABC.Score': abc_scores}
    df = pd.DataFrame(data)
    df['CellTypes'] = ['K562']+['other']*len(cells)
    print("Plotting the Pairplot")
    list_metrics = [dhs, h3k27ac, dhs_tss, h3k27ac_tss, promoter_activity]
    pearsonr, spearmanr= calculateMetrics(list_metrics)
    g = sns.pairplot(apply_jitter(df), x_vars={'logEnhancers.DHS.RPM', 'logEnhancers.H3K27ac.RPM', 'activity_base'}, y_vars={'logDHS.RPM.TSS1kb', 'logH3K27ac.RPM.TSS1kb', 'promoter_activity'}, hue='CellTypes', diag_kind='hist')
    g.savefig(os.path.join(files_outdir, "{}_{}.Correlation.png".format(name, target_gene)))
    g.savefig(os.path.join(outdir, "{}_{}.Correlation.pdf".format(name, target_gene)), format='pdf')
    plt.close('all')
    plt.close()

    with open(os.path.join(outdir, "{}_{}_celllines_correlation_values.txt".format(name, target_gene)), "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(pearsonr, spearmanr))
    print("Saving Plots...") 
    single_data = scale_values(df)
    single_data.to_csv(os.path.join(outdir, "{}_{}_ScaledForNull.tsv.gz".format(name, target_gene)), compression="gzip", sep="\t", index=False)
    ax = sns.pairplot(x_vars={'logEnhancers.DHS.RPM_median_centered', 'logEnhancers.H3K27ac.RPM_median_centered', 'activity_base_median_centered'}, y_vars={'logDHS.RPM.TSS1kb_median_centered', 'logH3K27ac.RPM.TSS1kb_median_centered', 'promoter_activity_median_centered'}, data=single_data, hue='CellTypes')
    ax.fig.set_size_inches(10,10)
    
    def plot_unity(xdata, ydata, **kwargs):
        mn = min(xdata.min(), ydata.min())
        mx = max(xdata.max(), ydata.max())
        points = np.linspace(mn, mx, 100)
        plt.gca().plot(points, points, color='k', marker=None,linestyle='--', linewidth=1.0)
    
    ax.map_offdiag(plot_unity)
    ax.savefig(os.path.join(files_outdir,"{}_{}_ScaledAndMedianCenteredPlot.png".format(name, target_gene)))
    plt.close("all")
    plt.close()
    
    df.to_csv(os.path.join(outdir, "{}_{}_celllines.tsv.gz".format(name, target_gene)), compression="gzip", sep="\t", index=False)
    
    return pearsonr_vals, spearmanr_vals, abc_scores, distance#, dist_not_found

