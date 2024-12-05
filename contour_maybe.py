# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 20:43:29 2024

@author: laure
"""
import seaborn as sns
import pandas
import numpy as np
import matplotlib.pyplot as plt

def true_series(parfile, num_exps):
  
  # create a df for the true values of the glitch
  if num_exps == 0:
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(parfile, sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']))
    return master_traits
  # recovery parameters if the par file is long enough
  if num_exps == 1:
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(parfile, sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_1"]['Value']))
    print("one exp. components found")
    return master_traits
                              
  if num_exps >= 2:
      cols = ["Element Name", "Value", "Fitting", "Error"]
      master_properties = pandas.read_csv(parfile, sep="\s+", names=cols)
      master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                      float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                      float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "GLF0D_1"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "GLTD_1"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "GLF0D_2"]['Value']),
                      float(master_properties.loc[master_properties['Element Name'] == "GLTD_2"]['Value']))
      
      print("two exp. components found")
      return master_traits
  
  else:
    print("no exp. components found")
    

    
def contour_plotter_5(data_name, no_of_exp, shade, master_traits, width,name):
    data = np.genfromtxt(data_name)
    
    sns.kdeplot(x=(data[:, 1] - master_traits[0] ), y= (data[:, 3] - master_traits[1]), levels=[0.01,0.33], ax = axs[0,0], color = shade, linestyles = width, alpha = [0.25,1])
    
    axs[0,0].set_xlabel(r"$\Delta \nu$ (Hz) - " + str(master_traits[1]), labelpad=15)
    axs[0,0].set_ylabel(r"$\Delta \dot \nu$ (Hz/s) - " + str(master_traits[3]), labelpad=15)
    
    
    if no_of_exp >= 1 :
        sns.kdeplot(x=(data[:, 11] - master_traits[6]), y=(data[:, 9] - master_traits[5]), levels=[0.01,0.33], ax = axs[0,1], color = shade, linestyles = width, alpha = [0.25,1])
        
        axs[0,1].set_xlabel(r"$\tau_{{d}}$("+str(1)+") - " + str(master_traits[6]) + "(days)", labelpad=15)
        axs[0,1].set_ylabel(r"$\Delta \nu_{{d}}$("+str(1)+") - " + str(master_traits[5]) + "(Hz)", labelpad=15)
        
        if no_of_exp == 2 :
            sns.kdeplot(x=(data[:, 19] - master_traits[8]), y=(data[:, 17] - master_traits[7]), levels=[0.01,0.33], ax = axs[0,2], color = shade, linestyles = width, alpha = [0.25,1],label = name)
            
            axs[0,2].set_xlabel(r"$\tau_{{d}}$("+str(2)+") - " + str(master_traits[8]) + "(days)", labelpad=15)
            axs[0,2].set_ylabel(r"$\Delta \nu_{{d}}$("+str(2)+") - " + str(master_traits[7]) + "(Hz)", labelpad=15)

def contour_plotter_15(data_name, no_of_exp, shade, master_traits,width):
    data = np.genfromtxt(data_name)
    
    sns.kdeplot(x=(data[:, 1] - master_traits[0] ), y= (data[:, 3] - master_traits[1]), levels=[0.01,0.33], ax = axs[1,0], color = shade, linestyles = width, alpha = [0.25,1])
    
    axs[1,0].set_xlabel(r"$\Delta \nu$ (Hz) - " + str(master_traits[1]), labelpad=15)
    axs[1,0].set_ylabel(r"$\Delta \dot \nu$ (Hz/s) - " + str(master_traits[3]), labelpad=15)
   
    
    if no_of_exp >= 1 :
        sns.kdeplot(x=(data[:, 11] - master_traits[6]), y=(data[:, 9] - master_traits[5]), levels=[0.01,0.33], ax = axs[1,1], color = shade, linestyles = width, alpha = [0.25,1])
        
        axs[1,1].set_xlabel(r"$\tau_{{d}}$("+str(1)+") - " + str(master_traits[6]) + "(days)", labelpad=15)
        axs[1,1].set_ylabel(r"$\Delta \nu_{{d}}$("+str(1)+") - " + str(master_traits[5]) + "(Hz)", labelpad=15)
        
        if no_of_exp == 2 :
            sns.kdeplot(x=(data[:, 19] - master_traits[8]), y=(data[:, 17] - master_traits[7]), levels=[0.01,0.33], ax = axs[1,2], color = shade, linestyles = width, alpha = [0.25,1])
            
            axs[1,2].set_xlabel(r"$\tau_{{d}}$("+str(2)+") - " + str(master_traits[8]) + "(days)", labelpad=15)
            axs[1,2].set_ylabel(r"$\Delta \nu_{{d}}$("+str(2)+") - " + str(master_traits[7]) + "(Hz)", labelpad=15)
            
def contour_plotter_30(data_name, no_of_exp, shade, master_traits,width):
    data = np.genfromtxt(data_name)
    
    sns.kdeplot(x=(data[:, 1] - master_traits[0] ), y= (data[:, 3] - master_traits[1]), levels=[0.01,0.33], ax = axs[2,0], color = shade, linestyles = width, alpha = [0.25,1])
    
    axs[2,0].set_xlabel(r"$\Delta \nu$ (Hz) - " + str(master_traits[1]), labelpad=15)
    axs[2,0].set_ylabel(r"$\Delta \dot \nu$ (Hz/s) - " + str(master_traits[3]), labelpad=15)
    
    if no_of_exp >= 1 :
        sns.kdeplot(x=(data[:, 11] - master_traits[6]), y=(data[:, 9] - master_traits[5]), levels=[0.01,0.33], ax = axs[2,1], color = shade, linestyles = width, alpha = [0.25,1])
        
        axs[2,1].set_xlabel(r"$\tau_{{d}}$("+str(1)+") - " + str(master_traits[6]) + "(days)", labelpad=15)
        axs[2,1].set_ylabel(r"$\Delta \nu_{{d}}$("+str(1)+") - " + str(master_traits[5]) + "(Hz)", labelpad=15)
        
        if no_of_exp == 2 :
            sns.kdeplot(x=(data[:, 19] - master_traits[8]), y=(data[:, 17] - master_traits[7]), levels=[0.01,0.33], ax = axs[2,2], color = shade, linestyles = width, alpha = [0.25,1])
            
            axs[2,2].set_xlabel(r"$\tau_{{d}}$("+str(2)+") - " + str(master_traits[8]) + "(days)", labelpad=15)
            axs[2,2].set_ylabel(r"$\Delta \nu_{{d}}$("+str(2)+") - " + str(master_traits[7]) + "(Hz)", labelpad=15)



glitch = 'c'

if glitch == 'b':

  master_par = "glitchB_master.par"

  arith_at_5 = "Glitch B @ 5\\arithmetic_1.5_glitchB_master_10000s_22_00.txt"
  arith_at_15 = "Glitch B @ 15\\arithmetic_4.33333_glitchB_master_10000s_13_42.txt"
  arith_at_30 = "Glitch B @ 30\logarithmic_35.2264_glitchB_master_10000s_14_15.txt"

  geo_at_5 = "Glitch B @ 5\geometric_1.6394_glitchB_master_10800s_18_25.txt"
  geo_at_15 = "Glitch B @ 15\geometric_1.66934_glitchB_master_10000s_15_38.txt"
  geo_at_30 = "Glitch B @ 30\logarithmic_35.2264_glitchB_master_10000s_14_15.txt"

  log_at_5 = "Glitch B @ 5\logarithmic_25.7197_glitchB_master_10500s_13_31.txt"
  log_at_15 = "Glitch B @ 15\logarithmic_34.76476_glitchB_master_10000s_12_17.txt"
  log_at_30 = "Glitch B @ 30\logarithmic_35.2264_glitchB_master_10000s_14_15.txt"

  peri_at_5 = "Glitch B @ 5\periodic_5_glitchB_master_10500s_13_31.txt"
  peri_at_15 = "Glitch B @ 15\periodic_15_glitchB_master_10350s_12_58.txt"
  peri_at_30 = "Glitch B @ 30\periodic_30_glitchB_master_10000s_14_15.txt"
  
  num_exps = 1

elif glitch == 'c':
  master_par = "glitch c @ 5\glitchC_master.par"

  arith_at_5 = "glitch c @ 5\\arithmetic_1.5_glitchC_master_10000s_12_41.txt"
  arith_at_15 = "glitch c @ 15\\arithmetic_4.33333_glitchC_master_10000s_16_04.txt"
  arith_at_30 = "glitch c @ 15\\arithmetic_4.33333_glitchC_master_10000s_16_04.txt"

  geo_at_5 = "glitch c @ 5\geometric_1.6394_glitchC_master_10800s_14_35.txt"
  geo_at_15 = "glitch c @ 15\geometric_1.6693_glitchC_master_10000s_18_15.txt"
  geo_at_30 = "glitch c @ 15\geometric_1.6693_glitchC_master_10000s_18_15.txt"
  
  log_at_5 = "glitch c @ 5\logarithmic_25.7197_glitchC_master_10500s_14_35.txt"
  log_at_15 = "glitch c @ 15\logarithmic_34.76476_glitchC_master_10000s_15_57.txt"
  log_at_30 = "glitch c @ 30\logarithmic_35.2264_glitchC_master_10000s_15_06.txt"
  
  peri_at_5 = "glitch c @ 5\periodic_5_glitchC_master_10500s_13_53.txt"
  peri_at_15 = "glitch c @ 15\periodic_15_glitchC_master_10350s_15_38.txt"
  peri_at_30 = "glitch c @ 30\periodic_30_glitchC_master_10000s_13_36.txt"

  num_exps = 2



#constants
counter = 0
counter_2 = 0
cadences = [5, 15, 30]
colours = ["tab:blue", "orange", "mediumorchid", "limegreen"]
linestyles = ["solid", "dashed", "dotted", "dashdot"]
name = ["arithmetic", "logarithmic", "geometric", "periodic"]

data_list_arith = [arith_at_5,arith_at_15,arith_at_30]
data_list_log = [log_at_5,log_at_15,log_at_30]
data_list_geo = [geo_at_5,geo_at_15,geo_at_30]
data_list_peri = [peri_at_5,peri_at_15,peri_at_30]

#creating figure
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(len(cadences), 1+num_exps, wspace = 0.5, hspace = 0.6)
axs = gs.subplots()

par = "glitch c @ 5\glitchC_master.par"

master = true_series(par,num_exps)

contour_plotter_5(arith_at_5,num_exps,colours[0],master,linestyles[0],name[0])
#contour_plotter_5(geo_at_5,num_exps,colours[1],master,linestyles[1],name[2])
#contour_plotter_5(peri_at_5,num_exps,colours[3],master,linestyles[2],name[3])
#contour_plotter_5(log_at_5,num_exps,colours[2],master,linestyles[3],name[1])

#contour_plotter_15(arith_at_15,num_exps,colours[0],master,linestyles[0])
#contour_plotter_15(geo_at_15,num_exps,colours[1],master,linestyles[1])
#contour_plotter_15(peri_at_15,num_exps,colours[3],master,linestyles[2])
#contour_plotter_15(log_at_15,num_exps,colours[2],master,linestyles[3])

#contour_plotter_30(arith_at_30,num_exps,colours[0],master,linestyles[0])
#contour_plotter_30(geo_at_30,num_exps,colours[1],master,linestyles[1])
#contour_plotter_30(peri_at_30,num_exps,colours[3],master,linestyles[2])
#contour_plotter_30(log_at_30,num_exps,colours[2],master,linestyles[3])

while counter <= 2:
    while counter_2 <= 2:
        axs[counter, counter_2].plot(0,0,marker = ".", color = "red")
        counter_2 = counter_2 + 1
    
    counter_2 = 0
    counter = counter + 1
    

ac5grp = fig.add_subplot(gs[0, 0:1])
ac5grp.set_yticks([])
ac5grp.set_xticks([])
ac5grp.set_frame_on(False)
ac15grp = fig.add_subplot(gs[1, 0:1])
ac15grp.set_yticks([])
ac15grp.set_xticks([])
ac15grp.set_frame_on(False)
ac30grp = fig.add_subplot(gs[2, 0:1])
ac30grp.set_yticks([])
ac30grp.set_xticks([])
ac30grp.set_frame_on(False)

sustaingrp = fig.add_subplot(gs[0:2, 0])
sustaingrp.set_yticks([])
sustaingrp.set_xticks([])
sustaingrp.set_frame_on(False)
recoverygrp = fig.add_subplot(gs[0:2, 1])
recoverygrp.set_yticks([])
recoverygrp.set_xticks([])
recoverygrp.set_frame_on(False)

#ac labels on the left  
ac5grp.set_ylabel("average cadenge: 5d", labelpad=50, size=10)
ac15grp.set_ylabel("average cadence: 15d", labelpad=50, size=10)
ac30grp.set_ylabel("average cadence: 30d", labelpad=50, size=10)

#sustain and recovery labels on the top
sustaingrp.tick_params(top = True, labeltop=True, labelbottom=False, bottom=False)
sustaingrp.set_title("Sustained Changes", pad=30, size=10)
recoverygrp.tick_params(labeltop=True, labelbottom=False, bottom=False, top = True)
recoverygrp.set_title("Recovery Components", pad=30, size=10)

#plt.legend(name,bbox_to_anchor=(1.05, 1.0), loc='upper left')

plt.show()
