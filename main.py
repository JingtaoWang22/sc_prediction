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












'''
GLOBAL VARIABLES:
'''
raw_dataset_folder = '../../cellar_data'



'''
HYPER PARAMETERS
'''





if __name__ == "__main__":

    # create dataset
    b,sc,genes,datanames  = create_dataset_folder(raw_dataset_folder)
