import pandas as pd
import numpy as np
import os

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def apply_jitter(data):
    for i in data.columns[:7]:
        data[i] = rand_jitter(data[i])
    return data

def make_dir(directory):
    if not os.path.exists(directory):
            os.makedirs(directory)
