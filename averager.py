# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:29:16 2024

@author: laure
"""

import numpy as np

def averager(data_name, num_of_exps):
    data = np.genfromtxt(data_name)
    
    GLF0_avg = np.average(data[:,1])
    GLF0_std = np.std(data[:,1])
    
    GLF1_avg = np.average(data[:,3])
    GLF1_std = np.std(data[:,3])
    
    if num_of_exps >= 1:
        if num_of_exps == 1:
            GLF0D_avg = np.average(data[:,9])
            GLF0D_std = np.std(data[:,9])
            
            GLTD_avg = np.average(data[:,11])
            GLTD_std = np.std(data[:,11])
            
            results = [GLF0_avg, GLF0_std, GLF1_avg, GLF1_std, GLF0D_avg, GLF0D_std, GLTD_avg, GLTD_std]
    
        if num_of_exps == 2 :
            GLF0D_avg = np.average(data[:,9])
            GLF0D_std = np.std(data[:,9])
            
            GLTD_avg = np.average(data[:,11])
            GLTD_std = np.std(data[:,11])
            
            GLF0D_2_avg = np.average(data[:,17])
            GLF0D_2_std = np.std(data[:,17])
            
            GLTD_2_avg = np.average(data[:,19])
            GLTD_2_std = np.std(data[:,19])
            
            results = [GLF0_avg, GLF0_std, GLF1_avg, GLF1_std, GLF0D_avg, GLF0D_std, GLTD_avg, GLTD_std, GLF0D_2_avg, GLF0D_2_std, GLTD_2_avg, GLTD_2_std]
            
        return results
            
glitch = 'b'

if glitch == 'b':

  master_par = "master files\glitchB_master.par"

  arith_at_5 = "Glitch B @ 5\\arithmetic_1.5_glitchB_master_10000s_22_00.txt"
  arith_at_15 = "Glitch B @ 15\\arithmetic_4.33333_glitchB_master_10000s_13_42.txt"
  arith_at_30 = "Glitch B @ 30\\arithmetic_1.8063_glitchB_master_10000s_16_49.txt"

  geo_at_5 = "Glitch B @ 5\geometric_1.6394_glitchB_master_10800s_18_25.txt"
  geo_at_15 = "Glitch B @ 15\geometric_1.66934_glitchB_master_10000s_15_38.txt"
  geo_at_30 = "Glitch B @ 30\geometric_3.5075_glitchB_master_10000s_16_51.txt"

  log_at_5 = "Glitch B @ 5\logarithmic_25.7197_glitchB_master_10500s_13_31.txt"
  log_at_15 = "Glitch B @ 15\logarithmic_34.76476_glitchB_master_10000s_12_17.txt"
  log_at_30 = "Glitch B @ 30\logarithmic_35.2264_glitchB_master_10000s_14_15.txt"

  peri_at_5 = "Glitch B @ 5\periodic_5_glitchB_master_10500s_13_31.txt"
  peri_at_15 = "Glitch B @ 15\periodic_15_glitchB_master_10350s_12_58.txt"
  peri_at_30 = "Glitch B @ 30\periodic_30_glitchB_master_10000s_15_10.txt"
  
  num_exps = 1
   
print(5)     
print(averager(peri_at_5,num_exps))
print(15)     
print(averager(peri_at_15,num_exps))
print(30)     
print(averager(peri_at_30,num_exps))