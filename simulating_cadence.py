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

    
    #print("calcing f0...", traits[0], master_traits[0])
    perc_f0 = (float(traits[0]) - master_traits[0])*100/master_traits[0] 
    #print(float(traits[0]))
    perc_f0_e = float(traits[1])/float(traits[0]) * perc_f0
    perc_f1 = (float(traits[2]) - master_traits[1])*100/master_traits[1] 
    perc_f1_e = float(traits[3])/float(traits[2]) * perc_f1

    ph = float(traits[4])
    
    print(perc_f0, perc_f0_e, perc_f1, perc_f1_e, ph)
    return perc_f0, perc_f0_e, perc_f1, perc_f1_e, ph

def tempo_nofit(par,tim):
    print("using tempo2 no fit")
    command_nofit = [
        "tempo2",
        "-f", par, tim,
        "-nofit",
        "noWarnings", ">&", "/dev/null",
        "-residuals"
        ]
    #print(' '.join(command_nofit), file=sys.stderr)
    proc = subprocess.Popen(command_nofit, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    out, err = proc.communicate()

def epoch_finder(par, tim):
    #runs tempo2 without a fit
    print("running tempo not fit")
    tempo_nofit(par, tim)
    #reads tempo2 generated residuals
    residuals = np.genfromtxt("residuals.dat")
    counter = 1
    error = 0.0001
    #finds estimation of glitch epoch
    while counter <= len(residuals):
        if np.abs(residuals[counter,1] - residuals[(counter -1),1]) > 10 * error:
            mid_point = (residuals[counter,0] + residuals[(counter -1),0])/2
            break 
            
        else :
            counter += 5
    
    return mid_point

def editting_par(parfile,GLEP,cols):
    new_line = np.empty(0)
    #reads in old par file
    lines = np.genfromtxt(parfile, delimiter="no-delim", dtype=str)
    for line in lines :
        if "GLEP_1" not in line :
            new_line = np.append(new_line,line)
      
    new_line = np.append(new_line,"GLEP_1          " + str(GLEP))    

    #saves it over the old par file
    np.savetxt(parfile, new_line, fmt="%s")
    

def run_fit(par, tim):
    command = [
        "tempo2",
        "-f", par, tim,
        "-nofit",
        "-fit", "GLF0_1",
        "-fit", "GLF1_1",
        "-fit", "GLPH_1",
        "-noWarnings",
        ]
    #print(' '.join(command), file=sys.stderr)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    out, err = proc.communicate()
    all_fields = out.split("\n")
    for this_field in all_fields:
        fields = this_field.split()
        print(fields)
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

def simulate(toas, sequence_type, const_args, sim_args, timfile, verbose = False, save_png = "results.png"):
    curr_iter = 0
    curr_sim_const = sim_args[0]
    step = np.abs(sim_args[1] - sim_args[0])/(sim_args[2]-1)
    steps = sim_args[2]
    master_par = "master_file.par"
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)
    master_traits = float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value'])

    print("running simulation for "+sequence_type+" sequence type\n[",end="")
    results = np.zeros((0,7))
    while curr_iter<steps:
        curr_iter += 1
        #print(toas)
        #print(indexes)
        # adds some 5d random variation so that we dont run into issues with the sample being the same every time
        start_randomiser = np.array([(const_args[1] + random.randint(0, 50)/10),
                            (const_args[1] + random.randint(0, 50)/10),
                            (const_args[1] + random.randint(0, 50)/10)])
        
        for offset in start_randomiser:
            print("inloop")
            passed_args = const_args[0], const_args[1]+offset, const_args[2], curr_sim_const
            indexes = tim_sampling.sample_from_toas(toas, sequence_type, passed_args, verbose)
            
            print("index array made")
            #if verbose: print(indexes)
            print(indexes)
            new_filename = "temp_toas.tim"
            print(toas)
            tim_sampling.gen_new_tim(timfile, indexes, new_filename)
            
            num_toas = len(indexes)
            #num_toas = tim_sampling.gen_new_tim(master_tim, indexes, new_filename)
            print("new toas generated, running tempo2")

            # run tempo2
            par, tim = "master_file_noglitch.par", new_filename
            
            min_MJD = round(np.min(toas))
            max_MJD = round(np.max(toas))
            
            initial_GLEP = random.randint(min_MJD,max_MJD)
            print(initial_GLEP)
            editting_par(par, initial_GLEP, cols)
            print("given par file initial guess")
            
            new_GLEP = epoch_finder(par, tim)
            print(new_GLEP)
            editting_par(par, new_GLEP, cols)
            print("given par accurate guess")
            
            print("running tempo2 with fit")
            traits = run_fit(par, tim)
            print(traits)
            
            print("retrieving results")
            print("master traits: ", master_traits)
            compare = compare_to_master(traits, master_traits)
            print(compare)
            curr_results = curr_sim_const, compare[0], compare[1], compare[2], compare[3], compare[4], num_toas
            print(curr_results)
            results = np.vstack((results, curr_results))
            
        print(str(curr_iter)+".", end="")
        sys.stdout.flush()
        #print("successfully simulated #"+ str(curr_iter)+ ", stepping log_const by "+str(step))
        curr_sim_const += step
        #print("the "+sequence_type+"_const is now "+str(curr_sim_const)+")")
    print("] - done!")
    # Below are settings used to generate a graphh.
    results = results.astype('float64')
    x = results[:,0].astype('float64')
    y = results[:,3].astype('float64')
    y_err = results[:,4]
    #z = np.log(results[:,4])
    #scaled_z = (z - z.min()) / z.ptp()
    #colours = plt.cm.Greys(scaled_z)
    
    plt.errorbar(x,np.abs(y),fmt=',')
    plt.scatter(x,np.abs(y),cmap='gist_gray',c=results[:,6],norm=colors.LogNorm(),edgecolors='gray')
    
    plt.colorbar(label="num. of ToAs")
    plt.xlabel(sequence_type+" constant")
    plt.ylabel("absolute value of % diff of retrieved and actual GLF0")
    #minimum = find_peaks(-np.abs(y), distance= 2000)
    #min_constant = minimum[0]
    #print(y[min_constant])
    #plt.scatter(x[min_constant],np.abs(y[min_constant]), marker="x", color = "red")
    plt.tight_layout()
    plt.savefig("figures/"+save_png, dpi=400)
    
    #plt.clf()
    #plt.plot(residuals[:,0], residuals[:,1])
    #plt.axvline(x = mid_point)
    #plt.savefig("results_22_10_24_2.png", dpi=400)
    return
    
def main():
    
    timfile = "master_toas.tim"
    toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])

    # User changeable 
    # 'logarithmic', 'arithmetic', 'geometric', 'periodic'
    SEQUENCE_TYPE = 'periodic'
    cadence_start = 0.5
    marker_offset = 0
    max_gap = 20
    verbose = False

    #simulation parameters
    num_iterations = 10

    ## LOGARITHMIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'logarithmic'
    log_const_min = 0
    log_const_max = 3

    ## ARITHMETIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'arithmetic'
    sequential_increase_min = 1.2 # num. of days per observation that time between observations increases by
    sequential_increase_max = 4 # num. of days per observation that time between observations increases by

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

    simulate(toas, SEQUENCE_TYPE, const_args, sim_args,timfile)

    #print("number of toas: " + str(len(indexes)))

if __name__ == "__main__":
    main()
        



