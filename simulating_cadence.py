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
import matplotlib.colors as colors
import sys
import random
from scipy.signal import find_peaks

def compare_to_master(traits, master_traits):
    # f0 % diff

    
    perc_f0 = (float(traits[0]) - master_traits[0])*100/master_traits[0] 
    perc_f0_e = float(traits[1]) * 100 / master_traits[0] 
    perc_f1 = (float(traits[2]) - master_traits[1])*100/master_traits[1] 
    perc_f1_e = float(traits[3]) * 100 / master_traits[0] 

    ph = float(traits[4])
    
    return perc_f0, perc_f0_e, perc_f1, perc_f1_e, ph

def run_fit(par, tim):
    command = [
        "tempo2",
        "-f", par, tim,
        "-nofit",
        "-fit", "GLF0_1",
        "-fit", "GLF1_1",
        "-fit", "GLPH_1",
        "-newpar", "-noWarnings", ">&", "/dev/null"
        ]
    print(' '.join(command), file=sys.stderr)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    out, err = proc.communicate()
    all_fields = out.split("\n")
    for this_field in all_fields:
        fields = this_field.split()
        #print(fields)
        if len(fields) > 2.0:
            if fields[0] == "GLF0_1":
                f0 = fields[2]
                f0_e = fields[3]
            if fields[0] == "GLF1_1":
                f1 = fields[2]
                f1_e = fields[3]
            if fields[0] == "GLPH_1":
                ph = fields[2]
    try:
        return f0, f0_e, f1, f1_e, ph
    except UnboundLocalError:
        return None

def simulate(toas, sequence_type, const_args, sim_args):
    curr_iter = 0
    curr_sim_const = sim_args[0]
    step = np.abs(sim_args[1] - sim_args[0])/(sim_args[2]-1)
    master_par = "master_file.par"
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)
    master_traits = float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value'])


    results = np.zeros((0,5))
    while curr_iter<sim_args[2]:
        curr_iter += 1
        #print(toas)
        #print(indexes)
        # adds some 5d random variation so that we dont run into issues with the sample being the same every time
        start_randomiser = np.array([(const_args[1] + random.randint(0, 50)/10),
                            (const_args[1] + random.randint(0, 50)/10),
                            (const_args[1] + random.randint(0, 50)/10)])
        
        for offset in start_randomiser:
            passed_args = const_args[0], const_args[1]+offset, const_args[2], curr_sim_const
            indexes = tim_sampling.sample_from_toas(toas, sequence_type, passed_args)
            
            print("index array made")   
            new_filename = "temp_toas.tim"
            
            num_toas = tim_sampling.gen_new_tim(timfile, indexes, new_filename)
            print("new toas generated, running tempo2")

            # run tempo2
            par, tim = "master_file_noglitch.par", new_filename
            traits = run_fit(par, tim)
            #print(traits)
            
            print("retrieving results")
            compare = compare_to_master(traits, master_traits)
            #print(compare)
            curr_results = curr_sim_const, compare[0], compare[1], compare[2], compare[3], compare[4], num_toas
            results = np.vstack((results, curr_results))
            
        print("successfully simulated #"+ str(curr_iter)+ ", stepping log_const by "+str(step))
        curr_sim_const += step
        print("the "+sequence_type+"_const is now "+str(curr_sim_const)+")")
    # Below are settings used to generate a graphh.
    results = results.astype('float64')
    #print(results)
    x = results[:,0].astype('float64')
    y = results[:,1].astype('float64')
    y_err = results[:,2]
    #z = np.log(results[:,4])
    #scaled_z = (z - z.min()) / z.ptp()
    #colours = plt.cm.Greys(scaled_z)
    plt.errorbar(x,np.abs(y),yerr=y_err,cmap='gist_gray',c=results[:,6],norm=colors.LogNorm(),edgecolors='gray')
    plt.colorbar(label="num. of ToAs")
    plt.xlabel(sequence_type+" constant")
    plt.ylabel("absolute value of % diff of retrieved and actual GLF0")
    minimum = find_peaks(-np.abs(y), distance= 2000)
    min_constant = minimum[0]
    print(y[min_constant])
    plt.scatter(x[min_constant],np.abs(y[min_constant]), marker="x", color = "red")
    plt.tight_layout()
    plt.savefig("results_15_20_24.png", dpi=400)
    return
    
timfile = "master_toas_2.tim"
toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])

# User changeable 
# 'logarithmic', 'arithmetic', 'geometric', 'periodic'
SEQUENCE_TYPE = 'periodic'
cadence_start = 0.5
marker_offset = 0
max_gap = 20
verbose = False

#simulation parameters
num_iterations = 50

## LOGARITHMIC - 
# these parameters are only used if SEQUENCE_TYPE is 'logarithmic'
log_const_min = 0
log_const_max = 3

## ARITHMETIC - 
# these parameters are only used if SEQUENCE_TYPE is 'arithmetic'
sequential_increase_min = 2 # num. of days per observation that time between observations increases by
sequential_increase_max = 2 # num. of days per observation that time between observations increases by

## GEOMETRIC - 
# these parameters are only used if SEQUENCE_TYPE is 'geometric'
multiplicative_increase_min = 1.2 # factor time between observations is multiplied by after an observation
multiplicative_increase_max = 4 # factor time between observations is multiplied by after an observation

## PERIODIC - 
# these parameters are only used if SEQUENCE_TYPE is 'periodic'
period_min = 0.5
period_max = 20


if SEQUENCE_TYPE == 'logarithmic':
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (log_const_min, log_const_max, num_iterations)
#    indexes = tim_sampling.sample_from_toas(toas, 'logarithmic', args, verbose)
elif SEQUENCE_TYPE == 'arithmetic':
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (sequential_increase_min, sequential_increase_max, num_iterations)
#    indexes = tim_sampling.sample_from_toas(toas, 'arithmetic', args, verbose)
elif SEQUENCE_TYPE == 'geometric':
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (multiplicative_increase_min, multiplicative_increase_max, num_iterations)
#    indexes = tim_sampling.sample_from_toas(toas, 'geometric', args, verbose)
elif SEQUENCE_TYPE == 'periodic':
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (period_min, period_max, num_iterations)
#    indexes = tim_sampling.sample_from_toas(toas, 'periodic', args, verbose)
else:
    print("invalid sequence type. doing nothing.")    

simulate(toas, SEQUENCE_TYPE, const_args, sim_args)

#print("number of toas: " + str(len(indexes)))

    



