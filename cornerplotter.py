# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:05:41 2024

@author: jade
"""

import numpy as np
import corner as cn


all_data = np.genfromtxt("700_tims_15_sims_18_41.txt", usecols=[1,3,8,13,15])

print(all_data)

cn.corner(all_data, labels=("glf0", "glf1", "epoch", "f0", "f1"))