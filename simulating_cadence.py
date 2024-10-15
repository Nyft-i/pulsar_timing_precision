# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:33:34 2024

@author: jade
"""

import tim_sampling
import numpy as np
import subprocess

timfile = "master_toas.tim"
toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])

# User changeable 
# 'logarithmic', 'arithmetic', 'geometric', 'periodic'
SEQUENCE_TYPE = 'periodic'
cadence_start = 0.5
marker_offset = 0
max_gap = 20
verbose = False

## LOGARITHMIC - 
# these parameters are only used if SEQUENCE_TYPE is 'logarithmic'
log_const = 0.4

## ARITHMETIC - 
# these parameters are only used if SEQUENCE_TYPE is 'arithmetic'
sequential_increase = 2 # num. of days per observation that time between observations increases by

## GEOMETRIC - 
# these parameters are only used if SEQUENCE_TYPE is 'geometric'
multiplicative_increase = 2 # factor time between observations is multiplied by after an observation

## PERIODIC - 
# these parameters are only used if SEQUENCE_TYPE is 'periodic'
period = 1

if SEQUENCE_TYPE == 'logarithmic':
    args = (cadence_start, marker_offset, max_gap, log_const)
    indexes = tim_sampling.sample_from_toas(toas, 'logarithmic', args, verbose)
elif SEQUENCE_TYPE == 'arithmetic':
    args = (cadence_start, marker_offset, max_gap, sequential_increase)
    indexes = tim_sampling.sample_from_toas(toas, 'arithmetic', args, verbose)
elif SEQUENCE_TYPE == 'geometric':
    args = (cadence_start, marker_offset, max_gap, multiplicative_increase)
    indexes = tim_sampling.sample_from_toas(toas, 'geometric', args, verbose)
elif SEQUENCE_TYPE == 'periodic':
    args = (cadence_start, marker_offset, max_gap, period)
    indexes = tim_sampling.sample_from_toas(toas, 'periodic', args, verbose)
else:
    print("invalid sequence type. doing nothing.")    


print("number of toas: " + str(len(indexes)))

new_filename = SEQUENCE_TYPE + "_toas.tim"
tim_sampling.gen_new_tim(timfile, indexes, new_filename)

par, tim = "master_file_noglitch.par", new_filename

subprocess.run([
    "tempo2",
    "-f", par, tim,
    "-nofit",
    "-fit", "GLF0_1",
    "-fit", "GLF1_1",
    "-fit", "GLPH_1",
    "-newpar"
    ])

properties = np.zeros((0,4))
properties = np.genfromtxt("new.par", skip_header=0, dtype=str)
print(properties)







