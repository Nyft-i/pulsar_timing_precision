# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:33:34 2024

@author: jade
"""

import tim_sampling
import numpy as np
import subprocess
import pandas

def compare_to_master(par, master_traits):
    # f0 % diff
    cols = ["Element Name", "Value", "Fitting", "Error"]
    properties = pandas.read_csv(par, sep="\s+", names=cols)
    
    perc_f0 = (properties.loc[properties['Element Name'] == "GLF0_1"]['Value'] - master_traits[0])*100/master_traits[0] 
    perc_f1 = (properties.loc[properties['Element Name'] == "GLF1_1"]['Value'] - master_traits[1])*100/master_traits[1] 
    ph = properties.loc[properties['Element Name'] == "GLPH_1"]['Value']
    
    return perc_f0, perc_f1, ph

    
    

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
master_par = "master_file.par"

subprocess.run([
    "tempo2",
    "-f", par, tim,
    "-nofit",
    "-fit", "GLF0_1",
    "-fit", "GLF1_1",
    "-fit", "GLPH_1",
    "-newpar"
    ])

cols = ["Element Name", "Value", "Fitting", "Error"]
master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)

master_traits = master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value'], master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value'], master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']

results = compare_to_master("new.par", master_traits)
print("results: ")
print(results[0])
print(results[1])
print(results[2])

    



