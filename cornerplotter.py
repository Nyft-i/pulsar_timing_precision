# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:05:41 2024

@author: jade
"""

import numpy as np
import corner as cn
import matplotlib.ticker as ticker
#from matplotlib import rcParams
import pandas

#"""
#rcParams["font.size"] = 16
#rcParams["font.family"] = "sans-serif"
#rcParams["font.sans-serif"] = ["Computer Modern Sans"]
#rcParams["text.usetex"] = True
#rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
#"""


all_data = np.genfromtxt("Glitch C @ 5/geometric_1.6394_glitchC_master_10800s_14_35.txt", usecols=[1,3,8,9,11,13,15,17,19])

cols = ["Element Name", "Value", "Fitting", "Error"]
master_properties = pandas.read_csv("glitchC_master.par", sep="\s+", names=cols)
master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "F0"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "F1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_2"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_2"]['Value']))

print(all_data[:,1]-master_traits[1])

all_data = np.array([
            all_data[:,0]-master_traits[0],
            all_data[:,1]-master_traits[1],
            all_data[:,2]-master_traits[4],
            all_data[:,3]-master_traits[5],
            all_data[:,4]-master_traits[6],
            all_data[:,5]-master_traits[7],
            all_data[:,6]-master_traits[8],
            all_data[:,7]-master_traits[9],
            all_data[:,8]-master_traits[10]
            ])

all_data = np.transpose(all_data)
            
print(all_data[1])

fig = cn.corner(all_data, labels=(r"$\Delta \nu \minus$"+str(master_traits[0]), 
                                  r"$\Delta \nu \minus$"+str(master_traits[1]),
                                  r"epoch$_g\minus$"+str(master_traits[3]), 
                                  r"$\Delta \nu_{{d,s}}\minus$"+str(master_traits[5]), 
                                  r"$\tau_{{d,s}}\minus$"+str(master_traits[6]), 
                                  r"$\nu \minus$"+str(master_traits[7]), 
                                  r"$\dot \nu \minus$"+str(master_traits[8]),
                                  r"$\Delta \nu_{{d,l}}\minus$"+str(master_traits[9]),
                                  r"$\tau_{{d,s}}\minus$"+str(master_traits[10])))


# get the axis
axes = np.array(fig.axes).reshape((9,9))
for i in range(9):
    for j in range(9):
      axes[i, j].set_xlabel(axes[i, j].get_xlabel(), fontsize=14)
      axes[i, j].set_ylabel(axes[i, j].get_ylabel(), fontsize=14)
      
fig.savefig("corner_plot.png", dpi=300)