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


    
def Param_Read(name,num_peaks=2):
    """opens .csv file of spectral parameters and unpacks data"""
    count = 0
    dictionary = {}
    p_dic = {}
    with open(name, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for row in reader:
            dictionary[str(count)] = row
            count += 1

    nums = []
    names = []
    for i in range(num_peaks):
      nums.append(int(0 + (5 * i)))
      names.append('peak_number_' + str(i))
      nums.append(int(1 + (5 * i)))
      names.append('peak_position_' + str(i))
      nums.append(int(2 + (5 * i)))
      names.append('peak_intensity_' + str(i))
      nums.append(int(3 + (5 * i)))
      names.append('peak_fwhm_' + str(i))
      nums.append(int(4 + (5 * i)))
      names.append('peak_area_' + str(i))
      
    nums.append(5 * num_peaks)
    names.append('index')
    for i in range(len(nums)):
        
        element_list = []
        for element in dictionary[str(nums[i])]:
            element_list.append(float(element))
            p_dic[names[i]] = element_list

    
  
    return p_dic

def value_correction(params,choice=1,upperbound='none',lowerbound='none'):
    """takes unpacked params, screens for bad data points and separates data
    into two groups based on the Kd classifier parameter"""
    
    
    
    
def mean(x):
    """computes the mean value of a given array"""
    
    elements = len(x)
    summ = sum(x)
    if elements == 0:
        a = 0
    else:
        a = summ/elements
    return a


    

def avg_params(name,params):
    """computes the average value for each spectral parameter"""
    
    dic1 = {}
    

    dic1['XD1'] = np.mean(params[0])
    dic1['sig0'] = np.std(params[0])
    dic1['WD1'] = np.mean(params[1])
    dic1['sig1'] = np.std(params[1])
    dic1['XG1'] = np.mean(params[2])
    dic1['sig2'] = np.std(params[2])
    dic1['WG1'] = np.mean(params[3])
    dic1['sig3'] = np.std(params[3])
    dic1['IDIG1'] = np.mean(params[4])
    dic1['sig4'] = np.std(params[4])
    dic1['KD1'] = np.mean(params[5])
    dic1['sig5'] = np.std(params[5])
    dic1['KG1'] = np.mean(params[6])
    dic1['sig6'] = np.std(params[6])
    dic1['DSLOPE1'] = np.mean(params[7])
    dic1['sig7'] = np.std(params[7])
    dic1['GSLOPE1'] = np.mean(params[8])
    dic1['sig8'] = np.std(params[8])
    dic1['T1'] = np.mean(params[9])
    dic1['sig9'] = np.std(params[9])
   
        
                      
     
param_Total = pd.DataFrame(param_T,index=names,columns=keyys)

#params_A.to_csv('A_TYPE_A_PARAMS.csv')
#params_B.to_csv('B_TYPE_B_PARAMS.csv')
#param_Total.to_csv('Total_PARAMS.csv')
