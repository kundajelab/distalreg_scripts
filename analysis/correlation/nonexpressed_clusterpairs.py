#! bin/python3 
from multiprocessing import Pool
import pandas as pd
import numpy as np 
import os, sys
import pyranges as pr
import pickle
import time
import gzip
import csv
import dask.dataframe as dd
from tqdm import tqdm

outdir = "/mnt/lab_data3/kmualim/PredictionFiles/"
save_outdir = "/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/fullk562_data/non_expressed"
cell_files = sys.argv[1]
cells = pd.read_csv(cell_files, sep="\t", header=None)

def process_cell(cell):
    print(cell)
    time.time()
    save_outdir = "/mnt/lab_data3/kmualim/CorrelationPlotsCellLines/splitK562_files/non_expressed"
    fpath = os.path.join(save_outdir, "{}_Intersected.tsv.gz".format(cell))
   
    columns={0:'chr', 1:'start', 2:'end', 3:'name', 4:'normalized_dhs', 5:'normalized_h3K27ac', 6:'activity_base', 7:'TargetGene', 8:'TargetGeneTSS', 9:'TargetGeneExpression', 10:'H3K27ac.RPM.TSS1Kb', 11:'DHS.RPM.TSS1Kb', 12:'distance', 13:'ABC.Score', 14:'chr_{}'.format(cell), 15:'start_{}'.format(cell), 16:'end_{}'.format(cell), 17:'name_{}'.format(cell), 19:'normalized_dhs_{}'.format(cell), 20:'normalized_h3K27ac_{}'.format(cell), 21:'activity_base_{}'.format(cell), 26:'TargetGene_{}'.format(cell), 27:'TargetGeneTSS_{}'.format(cell), 34:'Gene.H3K27ac.RPM.TSS1Kb', 36:'Gene.DHS.RPM.TSS1Kb', 38:'distance_{}'.format(cell), 47:'ABC.Score_{}'.format(cell)}

    total_cols = 52
    columns = ['num_{}'.format(str(idx)) if idx not in columns else columns[idx] for idx in range(total_cols)]
    data = pd.read_csv(fpath, sep="\t", compression="gzip", chunksize=10000000, names=columns)
    
    def include_row(row, i):
    #    if row['isSelfPromoter']:
    #        return False
    #    if row['TargetGene'].contains('LINCLOC'):
    #        return False
        if row['TargetGene'] != row['TargetGene_{}'.format(cell)]:
            return False
        if row['name']!=i:
            return False
        return True
    
    def keep_chunk(chunk):
        data_to_keep = chunk.loc[chunk['TargetGene']==chunk['TargetGene_{}'.format(cell)]]
        names = data_to_keep['name'].drop_duplicates()
        for name in names:
            outpath = '../non_expressed/files_v2/{}_Intersected.csv.gz'.format(name)
            chunkpername = data_to_keep.loc[data_to_keep['name']==name]
            chunkpername.to_csv(outpath, mode='a', sep="\t", compression="gzip")


    for chunk in tqdm(data):
        chunk['cellType'] = cell
        keep_chunk(chunk)
    for i in cells:
        yield i

if __name__=='__main__':
    generator = gen(cells[0])
    for cell in tqdm(generator):
        process_cell(cell)
    
    # To start preprocessing 
    #with Pool(15) as p:
    #    p.map(process_cell, list(cells[0]))

