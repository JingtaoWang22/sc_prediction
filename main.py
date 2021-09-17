#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 14:58:45 2021

@author: jingtao


handle training and testing the prediction
"""

import os
from utils import *
import numpy as np
import torch
from models import *

from sklearn.cluster import KMeans









'''
GLOBAL VARIABLES:
'''
raw_dataset_folder = '../../cellar_data'



'''
HYPER PARAMETERS
'''
PATIENT_BATCH = 1
N_PATIENT_CLUSTERS = 2




if __name__ == "__main__":

    # create dataset
    b,sc,genes,datanames  = create_dataset_folder(raw_dataset_folder)
    
    
    
    cluster_label = []
    modellist = []
    have_sc=[]
    '''
    for i in range(len(datanames)):
        model = initialize model
        modellist.append(model)
    '''
    
    # one iteration
    for iteration in range(2): # the iterations allowed by the budget given
        if (iteration == 0):
            #do clustering and testing in the first iteration
            kmeans = KMeans(n_clusters=N_PATIENT_CLUSTERS, random_state=0).fit(np.array(b))
            cluster_label = kmeans.labels_
            have_sc.append(0)
            
            
        
        for i in range(len(datanames)): # for each individual
            print("Trainign VAE for " + (i+1) + "th patient.")
            x = sc[i]
            bulk = b[i]
            genenames = genes[i]
        
            # fit the model if has sc

        ## choose the next batch of patients


# output simulations and errors



