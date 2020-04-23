#! bin/python3 

import pandas as pd
import numpy as np
import pickle
import argparse
import os
from joblib import Parallel, delayed, parallel_backend
from collections import ChainMap
from numpy import genfromtxt

def parse_args():
    parser = argparse.ArgumentParser(description="Making gene-centric dictionary")
##    parser.add_argument('--dataframe', help="Dataframe with scores")
    parser.add_argument('--weight', help="weight for promoter embeddings")
    parser.add_argument('--embeddings', help="Embeddings")
    parser.add_argument('--outdir', help="Outdir")
    parser.add_argument('--type', help="Type links")
    args = parser.parse_args()
    return args

# generator_function 
def gen_genes(genes):
    for gene in genes:
        yield gene

def generate_constant_promoter_weights(constantpromoterweight_df, weight):
    constantpromoterweight_df.iloc[:24672, 4] = weight
    return constantpromoterweight_df

def do_calc(dict_key, gene, jaspar_embeddings, distance):
    print("Gene: ", str(gene))
    indices = distance.loc[distance[3]==gene].index.astype(int)
    print(indices)
    total_dist_scores = [jaspar_embeddings[i] for i in indices]
    dict_key[gene] = np.sum(total_dist_scores, axis=0)
    return dict_key

def process(numpy_list_total, distance):
    dict_key = {}
    n_jobs=5
    genes = distance.iloc[:,3]
#    scores = np.array((distance.iloc[:,4]).ravel())
    with parallel_backend("loky", inner_max_num_threads=2):
        dict_key = Parallel(n_jobs=n_jobs)(delayed(do_calc)(dict_key, gene, numpy_list_total, distance) for gene in genes)
    return dict_key

def main(args, embeddings):
    print("Making directory")
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    print("Making constant promoter dataframe")
    constantpromoterweight_df = pd.read_csv("/mnt/lab_data2/kmualim/enh-gene-linking/datasets/embeddings/new_version/constant_promoter_weights/constantweightedpromoters_abc.csv", sep="\t", header=None)
    print("Adjusting weights for constant promoter weights dataframe")
    constantpromoterweight_df.iloc[:24672, 4] = args.weight
    print("Grabbing genes")
    genes = constantpromoterweight_df.iloc[:,3]
    print("Grabbing embeddings...")
    if embeddings.lower().find("basset") != -1: 
        embeddingfile = np.load(embeddings, mmap_mode='r')['embeddings']
    else:
        embeddingfile = np.load(embeddings, mmap_mode='r')['arr_0']
    print("Types of links")
    if args.type=='Random' or args.type=='all':
        print("Grabbing random links...")
        random_scores = np.random.rand(13995659, 1)
        calc_embeddings = [i*j for i,j in zip(random_scores, embeddingfile)]
        dict_key = process(calc_embeddings, constantpromoterweight_df)
        total_dict = dict(ChainMap(*dict_key))
        pickle.dump(total_dict, open(os.path.join(args.outdir, "RandomDictionary.p"), "wb"))
        del calc_embeddings
        del total_dict
    elif args.type=='Uniform' or args.type=='all':
        print("Grabbing uniform links...")
        calc_embedding = embeddingsfile
        dict_key = process(calc_embeddings, constantpromoterweight_df)
        total_dict = dict(ChainMap(*dict_key))
        pickle.dump(total_dict, open(os.path.join(args.outdir, "UniformDictionary.p"), "wb"))
        del calc_embeddings
        del total_dict
    elif args.type=='distance' or args.type=='all':
        print("Grabbing distance links...")
        df = pd.read_csv("/mnt/lab_data2/kmualim/enh-gene-linking/datasets/embeddings/new_version/constant_promoter_weights/Bassetdistanceweightedpromoters_abc.csv", sep="\t", header=None)
        random_scores = np.array((df.iloc[:,4]).ravel())
        calc_embeddings = [i*j for i,j in zip(random_scores, embeddingfile)]
        dict_key = process(calc_embeddings, df)
        total_dict = dict(ChainMap(*dict_key))
        pickle.dump(total_dict, open(os.path.join(args.outdir, "DistanceDictionary.p"), "wb"))
        del calc_embeddings
        del total_dict
    elif args.type=='abc' or args.type=='all':
        print("Grabbing abc links...")
        random_scores = np.array((constantpromoterweight_df.iloc[:,4]).ravel())
        calc_embeddings = [i*j for i,j in zip(random_scores, embeddingfile)]
        dict_key = process(calc_embeddings, df)
        total_dict = dict(ChainMap(*dict_key))
        pickle.dump(total_dict, open(os.path.join(args.outdir, "ABCDictionary.p"), "wb"))
        del calc_embeddings
        del total_dict
if __name__=="__main__":
    args = parse_args()
    main(args, args.embeddings)
