import pandas
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve

def true_series(parfile):
  num_exps = find_num_exps(parfile)
  
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
                              
  elif num_exps >= 2:
    master_traits = np.append(master_traits, (
                  float(par_read.loc[par_read['param'] == "GLF0D_2"]['value']),
                  float(par_read.loc[par_read['param'] == "GLTD_2"]['value'])
                  ))
    true_vals = pandas.DataFrame([master_traits], columns=["f0", "f1", "epoch", "phase", "df0", "df1", "df0d_1", "t_d_1", "df0d_2", "t_d_2"])
  
  

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
master_par = "glitchB_master.par"

arith_at_5 = "file"
arith_at_15 = "file"
arith_at_30 = "file"

geo_at_5 = "file"
geo_at_15 = "file"
geo_at_30 = "file"

log_at_5 = "file"
log_at_15 = "file"
log_at_30 = "file"

peri_at_5 = "file"
peri_at_15 = "file"
peri_at_30 = "file"

cols = ["const", "df0", "df0_e", "df1", "df1_e", "ph", "numtoas", "size", "epoch", "df0d_1", "df0d_e_1", "t_d_1", "t_d_e_1", "f0", "f0_e", "f1", "f1_e", "df0d_2", "df0d_2_e", "t_d_2", "t_d_2_e"]

data_arithmetic = np.array([pandas.read_csv(arith_at_5, names=cols, sep='\s+'), 
                            pandas.read_csv(arith_at_15, names=cols, sep='\s+'),
                            pandas.read_csv(arith_at_30, names=cols, sep='\s+')])

data_geometric = np.array([pandas.read_csv(geo_at_5, names=cols, sep='\s+'),
                            pandas.read_csv(geo_at_15, names=cols, sep='\s+'),
                            pandas.read_csv(geo_at_30, names=cols, sep='\s+')])

data_logarithmic = np.array([pandas.read_csv(log_at_5, names=cols, sep='\s+'),
                              pandas.read_csv(log_at_15, names=cols, sep='\s+'),
                              pandas.read_csv(log_at_30, names=cols, sep='\s+')])

data_periodic = np.array([pandas.read_csv(peri_at_5, names=cols, sep='\s+'), 
                          pandas.read_csv(peri_at_15, names=cols, sep='\s+'),
                          pandas.read_csv(peri_at_30, names=cols, sep='\s+')])

num_exps = find_num_exps(master_par)

# create a true values data series
true_vals = true_series(master_par)

cadences = [5, 15, 30]

# create a figure
fig, ax = plt.subplots(len(cadences), 1+num_exps, figsize=(8, 4))

for row, cadence in enumerate(cadences):
  print(cadence)
  
  # get the data for the current cadence
  data_arith = data_arithmetic[row]
  data_geo = data_geometric[row]
  data_log = data_logarithmic[row]
  data_peri = data_periodic[row]
  
  # create an array of the dataframes
  data_array = [data_arith, data_geo, data_log, data_peri]




  plt.sca(ax[row][0])

  colours = ["green", "purple", "black", "blue"]
  seq = ["Geometric", "Logarithmic", "Periodic"]

  # loop through the data array, df0 and df1
  for i, data in enumerate(data_array):  
    # create a 2d histogram from numpy to retrieve heights for a contour plot
    hist, xedges, yedges = np.histogram2d(data["df0"]-true_vals["df0"], data["df1"]-true_vals["df1"], bins=40)
    # numpy meshgrid
    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

    # astropy gaussian kernel smoothing (makes it not jagged and weird)
    kernel = Gaussian2DKernel(x_stddev=2)
    hist = convolve(hist, kernel)

    # plot the contour plot, the levels + contours surround the bins which contain more than that number of points i.e. where the density of points is higher than 10800/40^2
    cs = plt.contour(xmesh, ymesh, hist.T, levels=(1,10,30), colors=colours[i]) 
    plt.clabel(cs, inline=True, fontsize=10) # add labels to the contour lines
    
    # calculate the average and standard deviation of the data
    df0_avg = np.mean(data["df0"]-true_vals["df0"])
    df1_avg = np.mean(data["df1"]-true_vals["df1"])
    
    # plot the average and standard deviation
    plt.errorbar(df0_avg, df1_avg, fmt='x', color=colours[i], label=seq[i])
    
    # set the labels
    plt.xlabel(r"$\Delta \nu$ (Hz)")
    plt.ylabel(r"$\Delta \dot \nu$ (Hz/s)")
    plt.title("Sustained Changes")
    
  plt.scatter(0, 0, color="red", label="True Value") # plot the true value
    
    
  # loops for every recovery component
  for i in range(num_exps):
    plt.sca(ax[row][i+1])
    for j, data in enumerate(data_array):
      hist, xedges, yedges = np.histogram2d(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)], data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)], bins=40)
      xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])
      kernel = Gaussian2DKernel(x_stddev=2)
      hist = convolve(hist, kernel)
      cs = plt.contour(xmesh, ymesh, hist.T, levels=(1,10,30), colors=colours[j])
      plt.clabel(cs, inline=True, fontsize=10)
      plt.scatter(0, 0, color="red", label="True Value")
      
      #calculate the average and standard deviation of the data
      df0d_avg = np.mean(data["df0d_"+str(i+1)]-true_vals["df0d_"+str(i+1)])
      t_d_avg = np.mean(data["t_d_"+str(i+1)]-true_vals["t_d_"+str(i+1)])
      
      # plot the average and standard deviation
      plt.errorbar(t_d_avg, df0d_avg, fmt='x', color=colours[j], label=seq[j])
      
      # set the labels
      plt.xlabel(r"$\tau_{{d}}$("+str(i+1)+") (days)")
      plt.ylabel(r"$\Delta \nu_{{d}}$("+str(i+1)+") (Hz)")
      plt.title(f"Recovery Component {i+1}")

    
plt.show()

