#! bin/python3 

import pandas as pd
import numpy as np 
import os, sys
import pyranges as pr
import pickle
import time
import dask.dataframe as dd



def IntersectK562_CellTypes(cells, outdir):
    k562 = pd.read_csv("/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/fullk562_data/K562_EnhancerPredictionsAllPutative_noPromoterRegions.txt", sep="\t")
    k562_chr = k562.rename(columns={'chr':'Chromosome', 'start':'Start', 'end':'End'})
    k562_pyranges = pr.PyRanges(k562_chr)
    k562_clusters = k562_pyranges.cluster()    
    for i in cells[0]:
        t=time.time()
        print(i)
        print("Reading EnhancerPredictions...")
        data = pd.read_csv(os.path.join(outdir, "Predictions_{}_qnorm/EnhancerPredictionsAllPutative.txt.gz".format(i)), compression="gzip", sep="\t")
        data_ren = data.rename(columns={'chr':'Chromosome', 'start':'Start', 'end':'End', 'DHS.RPM':'Enhancers.DHS.RPM_{}'.format(i), 'distance':'distance_{}'.format(i),'H3K27ac.RPM':'Enhancers.H3K27ac.RPM_{}'.format(i), 'Gene.DHS.RPM.TSS1Kb':'DHS.RPM.TSS1Kb_{}'.format(i), 'Gene.H3K27ac.RPM.TSS1Kb':'H3K27ac.RPM.TSS1Kb_{}'.format(i), 'activity_base':'activity_base_{}'.format(i), 'TargetGenePromoterActivityQuantile':'TargetGenePromoterActivityQuantile_{}'.format(i), 'TargetGene':'TargetGene_{}'.format(i)})
        print("Makin pyranges")
        data_pyranges = pr.PyRanges(data_ren)
        print("Making pyranes cluster")
        data_clusters = data_pyranges.cluster()
        print("Intersecting K562 Pyranges clusters")
        intersected_dataframe = k562_clusters.join(data_clusters ,suffix="_{}".format(i))
        #subset_intersected_dataframe = intersected_dataframe.df
        #to_keep_df = subset_intersected_dataframe.loc[subset_intersected_dataframe['TargetGene'] == ['TargetGene_{}'.format(i)]]
        print("Saving file")
        intersected_dataframe.to_csv(os.path.join(save_outdir, "EnhancerPredictionsAllPutative_K562_{}_intersected.tsv.gz".format(i)), compression="gzip")
        print("Time taken: {}".format(str(time.time()-t)))

def process_cell(cell):
    time.time()
    j=cell
    save_outdir = "non_expressed"
    fpath = os.path.join(save_outdir, "{}_Intersected.tsv.gz".format(cell))
    data = pd.read_csv(fpath, compression="gzip", sep="\t", header=None)
    data_rem = data.rename(columns={0:'chr', 1:'start', 2:'end', 3:'name', 4:'normalized_dhs', 5:'normalized_h3K27ac', 6:'activity_base', 7:'TargetGene', 8:'TargetGeneTSS', 9:'TargetGeneExpression', 10:'H3K27ac.RPM.TSS1Kb', 11:'DHS.RPM.TSS1Kb', 12:'distance', 13:'ABC.Score', 14:'chr_{}'.format(j), 15:'start_{}'.format(j), 16:'end_{}'.format(j), 17:'name_{}'.format(j), 19:'normalized_dhs_{}'.format(j), 20:'normalized_h3K27ac_{}'.format(j), 21:'activity_base_{}'.format(j), 26:'TargetGene_{}'.format(j), 27:'TargetGeneTSS_{}'.format(j), 34:'Gene.H3K27ac.RPM.TSS1Kb', 36:'Gene.DHS.RPM.TSS1Kb', 38:'distance_{}'.format(j), 47:'ABC.Score_{}'.format(j) })
    intersected_data = data_rem.loc[data_rem['TargetGene']==data_rem['TargetGene_{}'.format(j)]]
    intersected_data_1 = intersected_data[['chr', 'start', 'end', 'name', 'normalized_dhs', 'normalized_h3K27ac', 'activity_base', 'TargetGene', 'TargetGeneExpression', 'H3K27ac.RPM.TSS1Kb','DHS.RPM.TSS1Kb','distance','ABC.Score', 'chr_{}'.format(j), 'start_{}'.format(j), 'end_{}'.format(j), 'name_{}'.format(j), 'normalized_dhs_{}'.format(j),'normalized_h3K27ac_{}'.format(j),'activity_base_{}'.format(j),'TargetGene_{}'.format(j),'TargetGeneTSS_{}'.format(j),'Gene.H3K27ac.RPM.TSS1Kb','Gene.DHS.RPM.TSS1Kb','distance_{}'.format(j),'ABC.Score_{}'.format(j)]]
    return intersected_data_1

#grab intersect for every cluster
def grab_indiv_cluster_intersect(file_a):
# open up cells 
    save_outdir="non_expressed"
    cells = pd.read_csv("SeventyOne_CellTypes.txt", sep="\t", header=None)
    dataframe_dict = pickle.load(open("Correlation_Intersected_ThresholdedEnhancers.p", "rb"))

    for celltype in cells[0]:
        print(celltype)
        t=time.time()

        print("Makin dataframe")
        t=time.time()
        dataframe_dict[celltype] = process_cell(celltype) 
        print("formed data")
        print("Time Taken:{}".format(str(time.time()-t)))
        print("Writing to HDF5 file")
    pickle.dump(dataframe_dict, open("Correlation_Intersected_ThresholdedEnhancers_v1.p", "wb"))
    print("Opened dictionary")
    
    clusters = pd.read_csv(file_a, sep="\t", header=None) 
    for name in clusters[0]:
        print(name)
        for j in cells[0]: 
            data = dataframe_dict[j]
            data = data.rename(columns={'Gene.H3K27ac.RPM.TSS1Kb':'Gene.H3K27ac.RPM.TSS1Kb_{}'.format(str(j)), 'Gene.DHS.RPM.TSS1Kb':'Gene.DHS.RPM.TSS1Kb_{}'.format(str(j))})
            # grab dataframe with similar clusters 
            matches = data.loc[data['name']==name]
            
            if j==cells[0][0]:
                df = matches
            else: 
                df = pd.concat([x, matches])
            x = df
        
        # fiure out GLS using this as well
        x.to_csv(os.path.join(save_outdir, "{}_Intersected.tsv.gz".format(name)), compression="gzip", sep="\t")

if __name__=='__main__':
    #IntersectK562_CellTypes(cells, outdir)
    file_a = sys.argv[1]
    grab_indiv_cluster_intersect(file_a)

