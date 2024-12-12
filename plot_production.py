import pandas
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
def true_series(parfile, num_exps):
  
  # create a df for the true values of the glitch
  par_read = pandas.read_csv = pandas.read_csv(parfile, names=["param", "value", "fitting", "error"], sep='\s+')
  master_traits = (float(par_read.loc[par_read['param'] == "F0"]['value']),
                  float(par_read.loc[par_read['param'] == "F1"]['value']),
                  float(par_read.loc[par_read['param'] == "GLEP_1"]['value']),
                  float(par_read.loc[par_read['param'] == "GLPH_1"]['value']),
                  float(par_read.loc[par_read['param'] == "GLF0_1"]['value']), 
                  float(par_read.loc[par_read['param'] == "GLF1_1"]['value']))
  true_vals = pandas.DataFrame([master_traits], columns=["f0", "f1", "epoch", "phase", "df0", "df1"])
  
  
  # recovery parameters if the par file is long enough
  if num_exps >= 1:
    master_traits = np.append(master_traits, (
                  float(par_read.loc[par_read['param'] == "GLF0D_1"]['value']),
                  float(par_read.loc[par_read['param'] == "GLTD_1"]['value'])
                  ))
    true_vals = pandas.DataFrame([master_traits], columns=["f0", "f1", "epoch", "phase", "df0", "df1", "df0d_1", "t_d_1"])
    print("one exp. components found")
                    
                              
  if num_exps >= 2:
    master_traits = np.append(master_traits, (
                  float(par_read.loc[par_read['param'] == "GLF0D_2"]['value']),
                  float(par_read.loc[par_read['param'] == "GLTD_2"]['value'])
                  ))
    true_vals = pandas.DataFrame([master_traits], columns=["f0", "f1", "epoch", "phase", "df0", "df1", "df0d_1", "t_d_1", "df0d_2", "t_d_2"])
    print("two exp. components found")
  
  else:
    print("no exp. components found")
    
  
  
  # create the dataframe
  true_vals = true_vals.squeeze(axis=0)
  return true_vals


def find_num_exps(parfile):
  # read the par file
  par_read = pandas.read_csv(parfile, names=["param", "value", "fitting", "error"], sep='\s+')
  if len(par_read) >= 17:
    return 2
  if len(par_read) >= 14:
    return 1
  else:
    return 0


# define columns for panda array
# according to this list of index titles:
"""
0 : const
1 : df0
2 : df0_e
3 : df1
4 : df1_e
5 : ph
6 : numtoas
7 : size
8 : epoch
9 : df0d_1
10: df0d_e_1
11: t_d_1
12: t_d_e_1
13: f0
14: f0_e
15: f1
16: f1_e
17: df0d_2
18: df0d_2_e
19: t_d_2
20: t_d_2_e
"""

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
  peri_at_30 = "Glitch B @ 30\logarithmic_35.2264_glitchB_master_10000s_14_15.txt"
  
  num_exps = 1
  

elif glitch == 'c':
  master_par = "glitch c\glitch c @ 5\glitchC_master.par"

  arith_at_5 = "glitch c\glitch c @ 5\\arithmetic_1.5_glitchC_master_10000s_12_41.txt"
  arith_at_15 = "glitch c\glitch c @ 15\\arithmetic_4.33333_glitchC_master_10000s_16_04.txt"
  arith_at_30 = "glitch c\glitch c @ 15\\arithmetic_4.33333_glitchC_master_10000s_16_04.txt"

  geo_at_5 = "glitch c\glitch c @ 5\geometric_1.6394_glitchC_master_10800s_14_35.txt"
  geo_at_15 = "glitch c\glitch c @ 15\geometric_1.6693_glitchC_master_10000s_18_15.txt"
  geo_at_30 = "glitch c\glitch c @ 15\geometric_1.6693_glitchC_master_10000s_18_15.txt"
  
  log_at_5 = "glitch c\glitch c @ 5\logarithmic_25.7197_glitchC_master_10500s_14_35.txt"
  log_at_15 = "glitch c\glitch c @ 15\logarithmic_34.76476_glitchC_master_10000s_15_57.txt"
  log_at_30 = "glitch c\glitch c @ 15\logarithmic_34.76476_glitchC_master_10000s_15_57.txt"
  
  peri_at_5 = "glitch c\glitch c @ 5\periodic_5_glitchC_master_10500s_13_53.txt"
  peri_at_15 = "glitch c\glitch c @ 15\periodic_15_glitchC_master_10350s_15_38.txt"
  peri_at_30 = "glitch c\glitch c @ 15\periodic_15_glitchC_master_10350s_15_38.txt"

  num_exps = 2


res2d = 20
contours= (0.02, 0.33)
cols = ["const", "df0", "df0_e", "df1", "df1_e", "ph", "numtoas", "size", "epoch", "df0d_1", "df0d_e_1", "t_d_1", "t_d_e_1", "f0", "f0_e", "f1", "f1_e", "df0d_2", "df0d_2_e", "t_d_2", "t_d_2_e"]

data_arithmetic = [pandas.read_csv(arith_at_5, names=cols, sep='\s+'), 
                            pandas.read_csv(arith_at_15, names=cols, sep='\s+'),
                            pandas.read_csv(arith_at_30, names=cols, sep='\s+')]

data_geometric = [pandas.read_csv(geo_at_5, names=cols, sep='\s+'),
                            pandas.read_csv(geo_at_15, names=cols, sep='\s+'),
                            pandas.read_csv(geo_at_30, names=cols, sep='\s+')]

data_logarithmic = [pandas.read_csv(log_at_5, names=cols, sep='\s+'),
                              pandas.read_csv(log_at_15, names=cols, sep='\s+'),
                              pandas.read_csv(log_at_30, names=cols, sep='\s+')]

data_periodic = [pandas.read_csv(peri_at_5, names=cols, sep='\s+'), 
                          pandas.read_csv(peri_at_15, names=cols, sep='\s+'),
                          pandas.read_csv(peri_at_30, names=cols, sep='\s+')]



# create a true values data series
true_vals = true_series(master_par, num_exps)

cadences = [5, 15, 30]

# create a figure
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(len(cadences), 1+num_exps, wspace = 0.4, hspace = 0.4)
ax = gs.subplots()


for row, cadence in enumerate(cadences):
  print("working at cadence", cadence)
  
  # get the data for the current cadence
  data_arith = data_arithmetic[row]
  data_geo = data_geometric[row]
  data_log = data_logarithmic[row]
  data_peri = data_periodic[row]
  
  # create an array of the dataframes
  data_array = [data_arith, data_geo, data_log, data_peri]




  plt.sca(ax[row][0])
  plt.scatter(0, 0, color="red", label="True Value", zorder = 5) # plot the true value


  colours = ["tab:blue", "orange", "mediumorchid", "limegreen"]
  seq = ["arithmetic", "geometric", "logarithmic", "periodic"]

  # loop through the data array, df0 and df1
  for i, data in enumerate(data_array):  
    # calculate the average and standard deviation of the data
    df0_avg = np.mean(data["df0"]-true_vals["df0"])
    df1_avg = np.mean(data["df1"]-true_vals["df1"])
    df0_std = np.std(data["df0"]-true_vals["df0"])
    df1_std = np.std(data["df1"]-true_vals["df1"])

    # removes all values which are more than 20 standard deviations away from the mean 
    badpoints = len(data) - len(data[(data["df0"]-true_vals["df0"] < 20*df0_std) & (data["df1"]-true_vals["df1"] > -20*df1_std)])
    if badpoints > 0: print(f"Removed {badpoints} points from {seq[i]} at {cadence} cadence")
    data = data[(data["df0"]-true_vals["df0"] < 20*df0_std) & (data["df1"]-true_vals["df1"] > -20*df1_std)]
    
    
    # create a 2d histogram from numpy to retrieve heights for a contour plot
    hist, xedges, yedges = np.histogram2d(data["df0"]-true_vals["df0"], data["df1"]-true_vals["df1"], bins=res2d)
    # numpy meshgrid
    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

    # astropy gaussian kernel smoothing (makes it not jagged and weird)
    kernel = Gaussian2DKernel(x_stddev=1)
    hist = convolve(hist, kernel)
    
    # rescale hist to be from 0 to 1
    hist = hist/np.max(hist)
    
    # plot the contour plot, the levels + contours surround the bins which contain more than that number of points i.e. where the density of points is higher than 10800/40^2
    cs = plt.contour(xmesh, ymesh, hist.T, levels=contours, colors=colours[i], alpha=0.5, zorder = 1) 
    fmt = {}
    for b in contours:
      fmt[b] = str(int((1-b)*100))+"%"
    #plt.clabel(cs, inline=True, fontsize=10, fmt=fmt) # add labels to the contour lines
    
    # plot the average and standard deviation
    plt.errorbar(df0_avg, df1_avg, fmt='x', color=colours[i], label=seq[i], zorder = 10)
    
    # plot each point with low alpha
    #plt.scatter(data["df0"]-true_vals["df0"], data["df1"]-true_vals["df1"], marker = ",", s = 1, alpha=0.1, color=colours[i], zorder = 0)
    
    
    # set the labels
    plt.xlabel(r"$\Delta \nu$ (Hz)")
    plt.ylabel(r"$\Delta \dot \nu$ (Hz/s)")
    
    
    
  # loops for every recovery component
  for i in range(num_exps):
    plt.sca(ax[row][i+1])
    plt.scatter(0, 0, color="red", label="True Value", zorder = 5) # plot the true value

    for j, data in enumerate(data_array):
      #calculate the average and standard deviation of the data
      df0d_avg = np.mean(data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)])
      t_d_avg = np.mean(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)])
      df0d_std = np.std(data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)])
      t_d_std = np.std(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)])
      
      # removes all values which are more than 20 standard deviations away from the mean
      data = data[(data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)] < 20*df0d_std) & (data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)] > -20*df0d_std)]
      
      # prints out the values where the timescale is close to 5
      print(data[(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)] < 0.1) & (data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)] > -0.1)])
      
      hist, xedges, yedges = np.histogram2d(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)], data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)], bins=res2d)
      xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])
      kernel = Gaussian2DKernel(x_stddev=1)
      hist = convolve(hist, kernel)
      
      # rescale hist to be from 0 to 1
      hist = hist/np.max(hist)
      
      # plotting
      cs = plt.contour(xmesh, ymesh, hist.T, levels=contours, colors=colours[j], alpha=0.5, zorder = 1)
      fmt = {}
      for b in contours:
        fmt[b] = str(int((1-b)*100))+"%"
      #plt.clabel(cs, inline=True, fontsize=10, fmt = fmt)
      plt.scatter(0, 0, color="red", label="True Value")
      
      # plot the average and standard deviation
      plt.errorbar(t_d_avg, df0d_avg, fmt='x', color=colours[j], label=seq[j], zorder = 10)
      
      # plot each point with low alpha
      #plt.scatter(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)], data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)], marker = ",", s = 1, alpha=0.1, color=colours[j], zorder = 0)
      
      # set the labels
      plt.xlabel(r"$\tau_{{d}}$("+str(i+1)+") (days)")
      plt.ylabel(r"$\Delta \nu_{{d}}$("+str(i+1)+") (Hz)")
      
      if i == num_exps-1 and row == len(cadences)-1:
        plt.legend()


# groups for labels

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


plt.show()

