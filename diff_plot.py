import numpy as np
import pandas
import tim_sampling
import simulating_cadence
import time

import matplotlib.pyplot as plt


# Plots our DDnu and DDnudot results for each of the cadence strategies
    
# Using pandas to read in the master file, probably a better way to do this but it works for now.
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
print(master_traits)


"""no recovery // gltichA"""
all_results_peri = np.genfromtxt("700_tims_15_sims_18_41.txt")
all_results_log = np.genfromtxt("140_tims_75_sims_18_55.txt")
all_results_geo = np.genfromtxt("90_tims_120_sims_13_14.txt")

"""recovery glitch b 5AC"""
#all_results_peri = np.genfromtxt("Glitch C @ 5/periodic_5_glitchC_master_10500s_13_53.txt")
#all_results_log = np.genfromtxt("Glitch C @ 5/logarithmic_25.7197_glitchC_master_10500s_14_35.txt")
#all_results_geo = np.genfromtxt("Glitch C @ 5/geometric_1.6394_glitchC_master_10800s_14_35.txt")

fig = plt.figure(figsize=(9, 3))
gs = fig.add_gridspec(1, 3, wspace = 0)
axs = gs.subplots(sharey = True, sharex = True)

fig.suptitle(r'difference in retrieved $\Delta \nu$ and $\Delta \dot \nu$ and actual values', x=0.5, y=1.05)
fig.supxlabel(r'$\Delta \nu - $' + str(master_traits[0]), y = -0.13)
fig.supylabel(r'$\Delta \dot \nu - $' + str(master_traits[1]), y=0.5, x=0.06)

# periodic
plt.sca(axs[0])
curr_data = all_results_peri
seq = "periodic"
col_a = "limegreen"
col_b = "darkgreen"
plt.scatter(curr_data[:,1]-master_traits[0], curr_data[:,3]-master_traits[1], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_f0 = np.mean(curr_data[:,1])
f0_std = np.std(curr_data[:,1])
avg_f1 = np.mean(curr_data[:,3])
f1_std = np.std(curr_data[:,3])
plt.errorbar(avg_f0-master_traits[0], avg_f1 - master_traits[1], xerr = f0_std, yerr = f1_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


# logarithmic
plt.sca(axs[1])
curr_data = all_results_log
seq = "logarithmic"
col_a = "mediumorchid"
col_b = "darkmagenta"
plt.scatter(curr_data[:,1]-master_traits[0], curr_data[:,3]-master_traits[1], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_f0 = np.mean(curr_data[:,1])
f0_std = np.std(curr_data[:,1])
avg_f1 = np.mean(curr_data[:,3])
f1_std = np.std(curr_data[:,3])
plt.errorbar(avg_f0-master_traits[0], avg_f1 - master_traits[1], xerr = f0_std, yerr = f1_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


# geometric
plt.sca(axs[2])
curr_data = all_results_geo
seq = "geometric"
col_a = "orange"
col_b = "darkorange"
plt.scatter(curr_data[:,1]-master_traits[0], curr_data[:,3]-master_traits[1], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_f0 = np.mean(curr_data[:,1])
f0_std = np.std(curr_data[:,1])
avg_f1 = np.mean(curr_data[:,3])
f1_std = np.std(curr_data[:,3])
plt.errorbar(avg_f0-master_traits[0], avg_f1 - master_traits[1], xerr = f0_std, yerr = f1_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)

axs[2].legend()

#set datetime in a format that can be used for saving the file
datetime = time.strftime("%Y-%m-%d-%H-%M")
plt.savefig("diff_plot-"+datetime+".png", dpi=400, bbox_inches="tight")
plt.show()

plt.clf()

fig = plt.figure(figsize=(9, 3))
gs = fig.add_gridspec(1, 3, wspace = 0)
axs = gs.subplots(sharey = True, sharex = True)

fig.suptitle(r'difference in retrieved recovery portion of $\Delta \nu$ and $\tau_r$ and actual values (short response)', x=0.5, y=1.05)
fig.supylabel(r'$\Delta \nu_{d,s}$ - ' + str(master_traits[5]), y=0.45, x=0.06)
fig.supxlabel(r'$\tau_{r,s}$ - ' + str(master_traits[6]), y = -0.05)

#periodic
plt.sca(axs[0])
curr_data = all_results_peri
seq = "periodic"
col_a = "limegreen"
col_b = "darkgreen"
plt.scatter(curr_data[:,11]-master_traits[6], curr_data[:,9]-master_traits[5], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,11])
x_std = np.std(curr_data[:,11])
avg_y = np.mean(curr_data[:,9])
y_std = np.std(curr_data[:,9])
plt.errorbar(avg_x-master_traits[6], avg_y - master_traits[5], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


#logarithmic
plt.sca(axs[1])
curr_data = all_results_log
seq = "logarithmic"
col_a = "mediumorchid"
col_b = "darkmagenta"
plt.scatter(curr_data[:,11]-master_traits[6], curr_data[:,9]-master_traits[5], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,11])
x_std = np.std(curr_data[:,11])
avg_y = np.mean(curr_data[:,9])
y_std = np.std(curr_data[:,9])
plt.errorbar(avg_x-master_traits[6], avg_y - master_traits[5], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


#geometric
plt.sca(axs[2])
curr_data = all_results_geo
seq = "geometric"
col_a = "orange"
col_b = "darkorange"
plt.scatter(curr_data[:,11]-master_traits[6], curr_data[:,9]-master_traits[5], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,11])
x_std = np.std(curr_data[:,11])
avg_y = np.mean(curr_data[:,9])
y_std = np.std(curr_data[:,9])
plt.errorbar(avg_x-master_traits[6], avg_y - master_traits[5], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)

#set datetime in a format that can be used for saving the file
datetime = time.strftime("%Y-%m-%d-%H-%M")
plt.savefig("diff_plot_rec_s-"+datetime+".png", dpi=400, bbox_inches="tight")
plt.show()

plt.clf()

fig = plt.figure(figsize=(9, 3))
gs = fig.add_gridspec(1, 3, wspace = 0)
axs = gs.subplots(sharey = True, sharex = True)

fig.suptitle(r'difference in retrieved recovery portion of $\Delta \nu$ and $\tau_r$ and actual values (long response)', x=0.5, y=1.05)
fig.supylabel(r'$\Delta \nu_{d,l}$ - ' + str(master_traits[9]), y=0.45, x=0.06)
fig.supxlabel(r'$\tau_{r,l}$ - ' + str(master_traits[10]), y = -0.05)


#periodic
plt.sca(axs[0])
curr_data = all_results_peri
seq = "periodic"
col_a = "limegreen"
col_b = "darkgreen"
plt.scatter(curr_data[:,19]-master_traits[10], curr_data[:,17]-master_traits[9], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,19])
x_std = np.std(curr_data[:,19])
avg_y = np.mean(curr_data[:,17])
y_std = np.std(curr_data[:,17])
plt.errorbar(avg_x-master_traits[10], avg_y - master_traits[9], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


#logarithmic
plt.sca(axs[1])
curr_data = all_results_log
seq = "logarithmic"
col_a = "mediumorchid"
col_b = "darkmagenta"
plt.scatter(curr_data[:,19]-master_traits[10], curr_data[:,17]-master_traits[9], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,19])
x_std = np.std(curr_data[:,19])
avg_y = np.mean(curr_data[:,17])
y_std = np.std(curr_data[:,17])
plt.errorbar(avg_x-master_traits[10], avg_y - master_traits[9], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)



#geometric
plt.sca(axs[2])
curr_data = all_results_geo
seq = "geometric"
col_a = "orange"
col_b = "darkorange"
plt.scatter(curr_data[:,19]-master_traits[10], curr_data[:,17]-master_traits[9], marker=',', zorder=1, alpha = 0.01, color = col_a)
avg_x = np.mean(curr_data[:,19])
x_std = np.std(curr_data[:,19])
avg_y = np.mean(curr_data[:,17])
y_std = np.std(curr_data[:,17])
plt.errorbar(avg_x-master_traits[10], avg_y - master_traits[9], xerr = x_std, yerr = y_std, zorder = 50, fmt = "x", color = col_b)
plt.title(seq)
plt.scatter(0, 0, c='r', label="real parameters", zorder =100)


#set datetime in a format that can be used for saving the file
datetime = time.strftime("%Y-%m-%d-%H-%M")
plt.savefig("diff_plot_rec_l-"+datetime+".png", dpi=400, bbox_inches="tight")
plt.show()


#fig.savefig("figures/fadbos.png", dpi=400, bbox_inches="tight")