# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:05:41 2024

@author: jade
"""

import numpy as np
import corner as cn
import matplotlib.ticker as ticker
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#from matplotlib import rcParams
import pandas

#"""
#rcParams["font.size"] = 16
#rcParams["font.family"] = "sans-serif"
#rcParams["font.sans-serif"] = ["Computer Modern Sans"]
#rcParams["text.usetex"] = True
#rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
#"""

""" glitch C
all_data = np.genfromtxt("Glitch C @ 5/geometric_1.6394_glitchC_master_10800s_14_35.txt", usecols=[1,3,8,9,11,13,15,17,19])
all_data = np.vstack((all_data, np.genfromtxt("Glitch C @ 5/arithmetic_1.5_glitchC_master_10000s_12_41.txt", usecols=[1,3,8,9,11,13,15,17,19])))
all_data = np.vstack((all_data, np.genfromtxt("Glitch C @ 5/periodic_5_glitchC_master_10500s_13_53.txt", usecols=[1,3,8,9,11,13,15,17,19])))
all_data = np.vstack((all_data, np.genfromtxt("Glitch C @ 5/logarithmic_25.7197_glitchC_master_10500s_14_35.txt", usecols=[1,3,8,9,11,13,15,17,19])))
"""


#all_data = np.genfromtxt("glitch c @ 15/arithmetic_4.33333_glitchC_master_10000s_16_04.txt", usecols=[1,3,8,9,11,13,15,17,19])

#glitch b data
all_data = np.genfromtxt("Glitch B @ 5/geometric_1.6394_glitchB_master_10800s_18_25.txt", usecols=[1,3,8,9,11,13,15])
all_data = np.vstack((all_data, np.genfromtxt("Glitch B @ 5/arithmetic_1.5_glitchB_master_10000s_22_00.txt", usecols=[1,3,8,9,11,13,15])))
all_data = np.vstack((all_data, np.genfromtxt("Glitch B @ 5/logarithmic_25.7197_glitchB_master_10500s_13_31.txt", usecols=[1,3,8,9,11,13,15])))
all_data = np.vstack((all_data, np.genfromtxt("Glitch B @ 5/periodic_5_glitchB_master_10500s_13_31.txt", usecols=[1,3,8,9,11,13,15])))



cols = ["Element Name", "Value", "Fitting", "Error"]
""" C
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
"""
# glitch B
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


print(all_data[:,1]-master_traits[1])

""" C
all_data = np.array([
            all_data[:,5]-master_traits[7],
            all_data[:,6]-master_traits[8],
            all_data[:,0]-master_traits[0],
            all_data[:,1]-master_traits[1],
            all_data[:,2]-master_traits[4],
            all_data[:,3]-master_traits[5],
            all_data[:,4]-master_traits[6],
            all_data[:,7]-master_traits[9],
            all_data[:,8]-master_traits[10]
            ])
"""

#GLITCH B
all_data = np.array([
            all_data[:,5]-master_traits[7],
            all_data[:,6]-master_traits[8],
            all_data[:,0]-master_traits[0],
            all_data[:,1]-master_traits[1],
            all_data[:,2]-master_traits[4],
            all_data[:,3]-master_traits[5],
            all_data[:,4]-master_traits[6]])


all_data = np.transpose(all_data)
            
print(all_data[1])

""" gltich C
fig = cn.corner(all_data, labels=(
                                  r"$\nu \minus$"+str(master_traits[7]), 
                                  r"$\dot \nu \minus$"+str(master_traits[8]),
                                  r"$\Delta \nu \minus$"+str(master_traits[0]), 
                                  r"$\Delta \dot \nu \minus$"+str(master_traits[1]),
                                  r"epoch$_g \minus$"+str(master_traits[3]), 
                                  r"$\Delta \nu_s \minus$"+str(master_traits[5]), 
                                  r"$\tau _s \minus$"+str(master_traits[6]), 
                                  r"$\Delta \nu_l \minus$"+str(master_traits[9]),
                                  r"$\tau_l \minus$ "+str(master_traits[10])), labelpad=0.1)
"""


#gltich B
fig = cn.corner(all_data, labels=(
                                  r"$\nu \minus$"+str(master_traits[7]), 
                                  r"$\dot \nu \minus$"+str(master_traits[8]),
                                  r"$\Delta \nu \minus$"+str(master_traits[0]), 
                                  r"$\Delta \dot \nu \minus$"+str(master_traits[1]),
                                  r"epoch$_g \minus$"+str(master_traits[3]), 
                                  r"$\Delta \nu_d \minus$"+str(master_traits[5]), 
                                  r"$\tau_d \minus$"+str(master_traits[6])), labelpad=0.1)


# get the axis
axes = np.array(fig.axes).reshape((7,7))
for i in range(7):
    for j in range(7):
      axes[i, j].set_xlabel(axes[i, j].get_xlabel(), fontsize=18)
      axes[i, j].set_ylabel(axes[i, j].get_ylabel(), fontsize=18)
      
#fig._gridspecs
#gs=fig.axes.get_gridspec()
#gs.update(hspace=0.1, wspace=0.1)

plt.subplots_adjust(wspace=0.16, hspace=0.16)
plt.show()
fig.savefig("corner_plot.png", dpi=400)