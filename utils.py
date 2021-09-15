#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:16:24 2021

@author: jingtao
"""

 
import anndata
import os
import torch 
import numpy as np
import pandas as pd
import timeit
import scanpy as sc

#datafolder = "newcastle74_normal"

folder = '../../cellar_data'
#dataname = "covid_portal_210320_with_raw.h5ad"

#adata = anndata.read("data/"+datafolder+'/'+dataname)

def load(path = '../../cellar_data/HBMP1_thymus_2_cellar.h5ad'):
    adata = anndata.read(path)
    return adata

def create_bulk(adata):
    '''
        crete bulk expression data according to sc data,
        assuming the matrix is (cell x gene).
    '''
    bulk = adata.X.mean(axis=0)
    adata.var['bulk'] = bulk
    return adata

def create_dataset_folder(folder_path=folder):
    '''
        create dataset given sc data folder
        
    '''
    if folder_path[-1] == '/':
        folder_path = folder_path[:-1]
        
    datalist = os.listdir(folder_path)
    b = []
    sc = []
    genes = []
    datanames = []
    for dataname in datalist:
        a = load(folder_path + '/' + dataname)
        if (type(a.X) != type(np.array([0]))):
            print(dataname+' skipped')
            # skip if no expression data or in sparse matrix format. will add support in future
            continue 
        a = create_bulk(a)
        b.append(np.array(a.var['bulk']))
        genes.append(list(a.var['bulk'].index))
        sc.append(a.X)
        datanames.append(dataname)
    return b,sc,genes,datanames

def preprocess(adata):
    satrt_time = timeit.default_timer()
    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    adata.var_names_make_unique()
    
    print(adata)
    
    #basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    #MT
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    sc.pp.scale(adata, max_value=10)
    
    end_time = timeit.defaulttimer()
    print("preprocessing took: "+str(end_time-start_time)+" seconds")
    
    return adata






if __name__ == "__main__":
    b,sc,genes=create_dataset_folder(folder)













