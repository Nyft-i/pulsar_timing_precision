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
import time

def compare_to_master(traits, master_traits):
    # f0 % diff

    
    #print("calcing f0...", traits[0], master_traits[0])
    perc_f0 = (float(traits[0]) - master_traits[0])*100/master_traits[0] 
    #print(float(traits))
    perc_f0_e = float(traits[1])/float(traits[0]) * perc_f0
    #print(perc_f0_e)
    perc_f1 = (float(traits[2]) - master_traits[1])*100/master_traits[1] 
    perc_f1_e = float(traits[3])/float(traits[2]) * perc_f1

    ph = float(traits[4])
    
    #print(perc_f0, perc_f0_e, perc_f1, perc_f1_e, ph)
    return perc_f0, perc_f0_e, perc_f1, perc_f1_e, ph, traits[5], traits[6]

def tempo_nofit(par,tim):
    #print("using tempo2 no fit")
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
    return np.genfromtxt("residuals.dat")

def epoch_finder(par, tim, master_traits):
    #runs tempo2 without a fit
    #print("running tempo not fit")
    residuals = tempo_nofit(par, tim)
    #reads tempo2 generated residuals
    counter = 1
    error = 0.0001
    #finds estimation of glitch epoch
    while counter <= len(residuals):
        if np.abs(residuals[counter,1] - residuals[(counter -1),1]) > 10 * error:
            change = ((residuals[counter,0] + residuals[(counter -1),0])/2) 
            mid_point = change + master_traits[3]
            break 
            
        else :
            counter += 5
    
    return mid_point

def editting_par(parfile,GLEP):
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
        "-noWarnings",">&","/dev/null"
        ]
    #print(' '.join(command), file=sys.stderr)
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
            if fields[0] == "MJD":
                epochs = fields[7], fields[9]
                epoch_e = fields[12]
    try:
        return f0, f0_e, f1, f1_e, ph, epochs, epoch_e
    except UnboundLocalError:
        return None

def simulate(toas, sequence_type, const_args, sim_args, verbose = False, master_tim="master_toas.tim", save_png = "results.png", num_sps=1):
    # This function samples TOAs from the master TOA file to a specific cadence strategy, then runs tempo2 on the new TOAs and compares the results to the master file.
    start_time = time.time()
    
    curr_iter = 0
    curr_sim_const = sim_args[0]
    step = np.abs(sim_args[1] - sim_args[0])/(sim_args[2]-1)
    steps = sim_args[2]
    master_par = "master_file.par"
    
    # Using pandas to read in the master file, probably a better way to do this but it works for now.
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']))

    print("running simulation for "+sequence_type+" sequence type\n[",end="")
    results = np.zeros((0,8))
    while curr_iter<steps:
        curr_iter += 1
        #print(toas)
        #print(indexes)
        
        # adds some 5d random variation so that we dont run into issues with the sample being the same every time
        start_randomiser = np.random.randint(0, 5, (num_sps))
            
        # For each offset, we generate a new set of toas, run tempo2, and compare the results to the master file
        for offset in start_randomiser:
            # We need passed args to take the form: cadence_start, offset, maxgap, const
            #const_args: start cadence, start offset, max_gap
            #sim_args: min const, max const, num_iterations
            
            passed_args = const_args[0], const_args[1]+offset, const_args[2], curr_sim_const
            
            print(passed_args)
            indexes = tim_sampling.sample_from_toas(toas, sequence_type, passed_args, verbose)
            num_toas = len(indexes)
            
            temp_tim = sequence_type+"_toas.tim"
            tim_sampling.gen_new_tim(master_tim, indexes, temp_tim)
            
            par, tim = "master_file_noglitch.par", temp_tim
            
            # Residual loading glep finder code, put it in the par file
            new_GLEP = epoch_finder(par, tim, master_traits)
            print(new_GLEP)
            editting_par(par, new_GLEP)
            
            # code for finding the closest TOA
            closest_index = (np.abs(master_traits[4] - toas[indexes]).argmin())
            distance_to_TOA = np.abs(master_traits[4] - toas[indexes][closest_index])
            #print(distance_to_TOA)
            
            # run tempo2
            traits = run_fit(par, tim)
            
            epochs = float(traits[5][0]), float(traits[5][1][:-1])
            closest_MJD_index = (np.abs(epochs - new_GLEP)).argmin()
            closest_MJD = epochs[closest_MJD_index]
            print(closest_MJD)
            
            # run tempo2 again with 0 phase MJD
            editting_par(par, closest_MJD)
            traits = run_fit(par, tim)
            
            # compare is an array of percentage differences between the retrieved and actual values of f0, f1, and phase (inc. error)
            compare = compare_to_master(traits, master_traits)
            curr_results = curr_sim_const, compare[0], compare[1], compare[2], compare[3], compare[4], num_toas, distance_to_TOA
            results = np.vstack((results, curr_results))
        
        # Print progress
        print(str(curr_iter)+"/"+str(steps)+" - curr " + str(sequence_type) + f" constant: {curr_sim_const:.2f}")
        sys.stdout.write("\033[K")
        sys.stdout.flush()
        curr_sim_const += step
        
    end_time = time.time()
    print("done! took " + f"{(end_time - start_time):.3f} seconds")
    
    # Below are settings used to generate a graphh.
    results = results.astype('float64')
    #print(results)
    x = results[:,0].astype('float64')
    y = results[:,1].astype('float64')
    y_err = results[:,2]
    
    plt.tight_layout()
    plt.errorbar(x,np.abs(y),xerr = 0, yerr = y_err,fmt=',')
    if sequence_type == 'logarithmic': ec = 'r'
    elif sequence_type == 'arithmetic': ec = 'b'
    elif sequence_type == 'geometric': ec = 'g'
    elif sequence_type == 'periodic': ec = 'y'
    plt.scatter(x,np.abs(y),cmap='gist_gray',c=results[:,6],s=results[:,7]*25,norm=colors.LogNorm(),edgecolors=ec)
    
    #plt.colorbar(label="num. of ToAs")
    #plt.xlabel(sequence_type+" constant")
    #plt.ylabel("absolute value of % diff of retrieved and actual GLF0_1")
    #plt.savefig("figures/glf0_"+save_png, dpi=400)
    #plt.clf()
    """
    y = results[:,3].astype('float64')
    y_err = results[:,4]

    plt.tight_layout()
    plt.errorbar(x,np.abs(y),xerr = 0, yerr = y_err,fmt=',')
    plt.scatter(x,np.abs(y),cmap='gist_gray',c=results[:,6],s=results[:,7]*25,norm=colors.LogNorm(),edgecolors='gray')
    
    plt.colorbar(label="num. of ToAs")
    plt.xlabel(sequence_type+" constant")
    plt.ylabel("absolute value of % diff of retrieved and actual GLF1_1")
    plt.savefig("figures/glf1_"+save_png, dpi=400)
    """
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
    num_iterations = 20

    ## LOGARITHMIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'logarithmic'
    log_const_min = 0.5
    log_const_max = 2

    ## ARITHMETIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'arithmetic'
    sequential_increase_min = 0.5 # num. of days per observation that time between observations increases by
    sequential_increase_max = 4 # num. of days per observation that time between observations increases by

    ## GEOMETRIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'geometric'
    multiplicative_increase_min = 1 # factor time between observations is multiplied by after an observation
    multiplicative_increase_max = 4 # factor time between observations is multiplied by after an observation

    ## PERIODIC - 
    # these parameters are only used if SEQUENCE_TYPE is 'periodic'
    period_min = 0.5
    period_max = 20


    SEQUENCE_TYPE = 'logarithmic'
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (log_const_min, log_const_max, num_iterations)
    simulate(toas, SEQUENCE_TYPE, const_args, sim_args,timfile)
    #    indexes = tim_sampling.sample_from_toas(toas, 'logarithmic', args, verbose)
    SEQUENCE_TYPE = 'arithmetic'
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (sequential_increase_min, sequential_increase_max, num_iterations)
    simulate(toas, SEQUENCE_TYPE, const_args, sim_args,timfile)
#    indexes = tim_sampling.sample_from_toas(toas, 'arithmetic', args, verbose)
    SEQUENCE_TYPE = 'geometric'
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (multiplicative_increase_min, multiplicative_increase_max, num_iterations)
    simulate(toas, SEQUENCE_TYPE, const_args, sim_args,timfile)
#    indexes = tim_sampling.sample_from_toas(toas, 'geometric', args, verbose)
    SEQUENCE_TYPE = 'periodic'
    const_args = (cadence_start, marker_offset, max_gap)
    sim_args = (period_min, period_max, num_iterations)
    simulate(toas, SEQUENCE_TYPE, const_args, sim_args,timfile)
#   indexes = tim_sampling.sample_from_toas(toas, 'periodic', args, verbose)
    
    plt.colorbar(label="num. of ToAs")
    plt.xlabel("sequence constant")
    plt.ylabel("absolute value of % diff of retrieved and actual GLF0_1")
    


    #print("number of toas: " + str(len(indexes)))

if __name__ == "__main__":
    main()
        



