# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 13:13:24 2019

@author: Jordan
"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
import pandas as pd

def num_name_gen(num_peaks):
    """generates names for your parameter dictionary based 
       on the number of peaks you're interested in"""
    
    nums = []
    names = []
    
    for i in range(num_peaks):
        nums.append(int(0 + (5 * i)))
        names.append('peak_number_' + str(i+1))
        nums.append(int(1 + (5 * i)))
        names.append('peak_position_' + str(i+1))
        nums.append(int(2 + (5 * i)))
        names.append('peak_intensity_' + str(i+1))
        nums.append(int(3 + (5 * i)))
        names.append('peak_fwhm_' + str(i+1))
        nums.append(int(4 + (5 * i)))
        names.append('peak_area_' + str(i+1))
    
    return names,nums

def Param_Read(file_name,num_peaks=2):
    """opens .csv file of spectral parameters and unpacks data"""
    
    if '.csv' in name:
        continue
    else:
        file_name = file_name + '.csv' 
        
    count = 0
    dictionary = {}
    params = {}
    
    
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for row in reader:
            dictionary[str(count)] = row
            count += 1
    
    names,nums = num_name_gen(num_peaks)
   
      
    nums.append(5 * num_peaks)
    names.append('index')
    for i in range(len(nums)):
        
        element_list = []
        for element in dictionary[str(nums[i])]:
            element_list.append(float(element))
            params[names[i]] = element_list

    
  
    return params

   
    
   
def avg_params(params):
    """computes the average value for each spectral parameter"""
    avg_spec_params = {}
    for key, value in params.items()
        avg_spec_params[key] = np.mean(np.array(value))
   
    return avg_spec_params
                      
     


