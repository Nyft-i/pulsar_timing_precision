# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:33:34 2024

@author: jade
"""

import tim_sampling
import numpy as np
import subprocess
import pandas
import matplotlib.pyplot as plt

def compare_to_master(par, master_traits):
    # f0 % diff
    cols = ["Element Name", "Value", "Fitting", "Error"]
    properties = pandas.read_csv(par, sep="\s+", names=cols)
    
    perc_f0 = (float(properties.loc[properties['Element Name'] == "GLF0_1"]['Value']) - master_traits[0])*100/master_traits[0] 
    perc_f1 = (float(properties.loc[properties['Element Name'] == "GLF1_1"]['Value']) - master_traits[1])*100/master_traits[1] 
    ph = float(properties.loc[properties['Element Name'] == "GLPH_1"]['Value'])
    
    return perc_f0, perc_f1, ph

def simulate(toas, sequence_type, const_args, sim_args):
    curr_iter = 0
    curr_log_const = sim_args[0]
    step = np.abs(sim_args[1] - sim_args[0])/(sim_args[2]-1)
    master_par = "master_file.par"
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)
    master_traits = float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value'])

    
    if sequence_type == 'logarithmic':
        results = np.zeros((0,4))
        while curr_iter<sim_args[2]:
            curr_iter += 1
            passed_args = const_args[0], const_args[1], const_args[2], curr_log_const
            #print(toas)
            indexes = tim_sampling.sample_from_toas(toas, sequence_type, passed_args)
            print("index array made")
            #print(indexes)
            new_filename = sequence_type + "_toas.tim"
            tim_sampling.gen_new_tim(timfile, indexes, new_filename)
            print("new toas generated, running tempo2")

            # run tempo2
            par, tim = "master_file_noglitch.par", new_filename

            subprocess.run([
                "tempo2",
                "-f", par, tim,
                "-nofit",
                "-fit", "GLF0_1",
                "-fit", "GLF1_1",
                "-fit", "GLPH_1",
                "-newpar", "-noWarnings"
                ])

            print("retrieving results")
            compare = compare_to_master("new.par", master_traits)
            curr_results = curr_log_const, compare[0], compare[1], compare[2]
            results = np.vstack((results, curr_results))
            print("successfully simulated #"+ str(curr_iter)+ ", stepping log_const by "+str(step))
            curr_log_const += step
            print("(log_const is now "+str(curr_log_const)+")")
        print(results)
        plt.plot(results[:,0], results[:,1])
        plt.savefig("results 15/20/24.png", dpi=400)
        return
    else:
        print("you should never see this")
    
timfile = "master_toas.tim"
toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])

# User changeable 
# 'logarithmic', 'arithmetic', 'geometric', 'periodic'
SEQUENCE_TYPE = 'logarithmic'
cadence_start = 0.5
marker_offset = 0
max_gap = 20
verbose = False

#simulation parameters
log_const_min = 0
log_const_max = 1
num_iterations = 100

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
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (log_const_min, log_const_max, num_iterations)
#    indexes = tim_sampling.sample_from_toas(toas, 'logarithmic', args, verbose)
elif SEQUENCE_TYPE == 'arithmetic':
    args = (cadence_start, marker_offset, max_gap, sequential_increase)
 #   indexes = tim_sampling.sample_from_toas(toas, 'arithmetic', args, verbose)
elif SEQUENCE_TYPE == 'geometric':
    args = (cadence_start, marker_offset, max_gap, multiplicative_increase)
  #  indexes = tim_sampling.sample_from_toas(toas, 'geometric', args, verbose)
elif SEQUENCE_TYPE == 'periodic':
    args = (cadence_start, marker_offset, max_gap, period)
   # indexes = tim_sampling.sample_from_toas(toas, 'periodic', args, verbose)
else:
    print("invalid sequence type. doing nothing.")    

simulate(toas, SEQUENCE_TYPE, const_args, sim_args)

#print("number of toas: " + str(len(indexes)))

    



