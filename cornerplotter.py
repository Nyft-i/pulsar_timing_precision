# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:05:41 2024

@author: jade
"""

import numpy as np
import corner as cn
import matplotlib.ticker as ticker


all_data = np.genfromtxt("logarithmic_25.7197_glitchB_master_10500s_13_31.txt", usecols=[1,3,8,10,12,13,15])

print(all_data)

fig = cn.corner(all_data, labels=("glf0", "glf1", "epoch", "glf0d", "gltd", "f0", "f1"))

ax = fig.axes
print(ax.shape())
for curr_ax in ax:
    curr_ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    curr_ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))