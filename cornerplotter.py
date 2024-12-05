# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:05:41 2024

@author: jade
"""

import numpy as np
import corner as cn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas


all_data = np.genfromtxt("geometric_1.6394_glitchB_master_10800s_18_25.txt", usecols=[1,3,8,9,11,13,15])

cols = ["Element Name", "Value", "Fitting", "Error"]
master_properties = pandas.read_csv("glitchB_master.par", sep="\s+", names=cols)
master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "F0"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "F1"]['Value']))
                    #float(master_properties.loc[master_properties['Element Name'] == "GLF0D_2"]['Value']),
                    #float(master_properties.loc[master_properties['Element Name'] == "GLTD_2"]['Value']))

print(all_data[:,1]-master_traits[1])

all_data = np.array([all_data[:,0]-master_traits[0],
            all_data[:,1]-master_traits[1],
            all_data[:,2]-master_traits[4],
            all_data[:,3]-master_traits[5],
            all_data[:,4]-master_traits[6],
            all_data[:,5]-master_traits[7],
            all_data[:,6]-master_traits[8]])
            #all_data[:,7]-master_traits[9],
            #all_data[:,8]-master_traits[10]])

all_data = np.transpose(all_data)
            
print(all_data[1])

plt.rcParams.update({'font.size': 11})

fig = cn.corner(all_data, labels=(f"glf0$\minus$"+str(master_traits[0]), 
                                  f"glf1$\minus$"+str(master_traits[1]),
                                  f"epoch$\minus$"+str(master_traits[3]), 
                                  f"glf0d_1$\minus$"+str(master_traits[5]), 
                                  f"gltd_1$\minus$"+str(master_traits[6]), 
                                  f"f0$\minus$"+str(master_traits[5]), 
                                  f"f1$\minus$"+str(master_traits[6])))
                                  #f"glf0d_2$\minus$"+str(master_traits[9]),
                                  #f"gltd_2$\minus$"+str(master_traits[10])))


ax = fig.axes