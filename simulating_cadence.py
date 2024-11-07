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
    # sort the residuals by the mjd
    residuals = residuals[residuals[:,0].argsort()]
    print(residuals[:,0])
    #reads tempo2 generated residuals
    counter = 1
    error = 0.0001
    #finds estimation of glitch epoch
    while counter <= len(residuals):
        difference = np.abs(residuals[counter,1] - residuals[(counter -1),1])
        #print("hello")
        #print(counter, difference)
        if difference > 30 * error:
            #print(residuals[counter,0])
            #print(residuals[counter-1,0])
            change = ((residuals[counter,0] + residuals[(counter -1),0])/2) 
            mid_point = change + master_traits[3]
            #plt.scatter(residuals[:,0], residuals[:,1])
            #plt.axvline(change, color = "pink")
            #plt.savefig("figures/midpoint_checker.png", dpi=400, bbox_inches="tight")
            #plt.clf()
            print(mid_point)
            break 
            
        else :
            counter += 1
    
    return mid_point

def editting_par(parfile, param, g_property="GLEP_1"):
    new_line = np.empty(0)
    #reads in old par file
    lines = np.genfromtxt(parfile, delimiter="no-delim", dtype=str)
    for line in lines :
        if g_property not in line :
            new_line = np.append(new_line, line) 
      
    new_line = np.append(new_line,g_property + " " + str(param))    

    #saves it over the old par file
    np.savetxt(parfile, new_line, fmt="%s")    

def run_fit(par, tim, recovery_mode = False, no_phase_fit = False):
    command = [
        "tempo2",
        "-f", par, tim,
        "-nofit",
        "-fit", "GLF0_1",
        "-fit", "GLF1_1",
        "-fit", "GLPH_1",
        "-fit", "F1",
        "-fit", "F0",
        "-noWarnings",">&","/dev/null"
        ]
    
    if no_phase_fit == True:
        command = [
            "tempo2",
            "-f", par, tim,
            "-nofit",
            "-fit", "GLF0_1",
            "-fit", "GLF1_1",
            "-fit", "F1",
            "-fit", "F0"
            ]
    
    if recovery_mode == True :
        command_rec = ["-fit", "GLF0D_1",
                       "-fit", "GLTD_1"]
        command = np.hstack((command,command_rec))
        #print(command)
    
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
            if fields[0] == "GLF0D_1":
                recovered_F0 = fields[2]
                recovered_F0_e = fields[3]
            if fields[0] == "GLTD_1":
                recovered_timescale = fields[2]
                recovered_timescale_e = fields[3]

            
    try:
        if recovery_mode == True:
            return f0, f0_e, f1, f1_e, ph, epochs, epoch_e, recovered_F0, recovered_F0_e, recovered_timescale, recovered_timescale_e
        
        return f0, f0_e, f1, f1_e, ph, epochs, epoch_e
    except UnboundLocalError:
        return None
    
def single_simulate(toas, sequence_type, const_args, sim_arg, verbose = False, master_tim="master_toas.tim", master_par="master_file.par", num_sps=1, epoch_finding_mode=False):
    # This function samples TOAs from the master TOA file to a specific cadence strategy, then runs tempo2 on the new TOAs and compares the results to the master file.
    start_time = time.time()
    
    # Using pandas to read in the master file, probably a better way to do this but it works for now.
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(master_par, sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']))
    
    # adds some 5d random variation so that we dont run into issues with the sample being the same every time
    passed_args = const_args[0], const_args[1], const_args[2], sim_arg
    strategy_period, strat_toas = tim_sampling.find_sequence_period_info(sequence_type, passed_args)
    start_randomiser = np.random.randint(0, strategy_period*10, (num_sps))
    start_randomiser = start_randomiser/10
    all_results = np.zeros((0,9))
    all_epochs = np.zeros(0)
    
    # For each offset, we generate a new set of toas, run tempo2, and compare the results to the master file
    print("running simulation for "+sequence_type+" sequence type\n[",end="")
    for number, offset in enumerate(start_randomiser):
        # We need passed args to take the form: cadence_start, offset, maxgap, const
        # const_args: start cadence, start offset, max_gap
        print("offset: ", offset, end=" ")
        passed_args = const_args[0], const_args[1]+offset, const_args[2], sim_arg
        
        indexes, num_toas = tim_sampling.sample_from_toas(toas, sequence_type, passed_args, verbose, strat_period=strategy_period)
        
        temp_tim = sequence_type+"_toas.tim"
        #print(indexes)
        tim_sampling.gen_new_tim(master_tim, indexes, temp_tim)
        
        par, tim = "master_noglitch_exp.par", temp_tim
        
        # Residual loading glep finder code, put it in the par file
        new_GLEP = epoch_finder(par, tim, master_traits)
        print(new_GLEP)
        editting_par(par, new_GLEP)
        
        # run tempo2
        traits = run_fit(par, tim)
        print("post tempo2 1 traits")
        print(traits)
        editting_par(par, traits[0], "GLF0_1")
        editting_par(par, traits[2], "GLF1_1")
        
        epochs = float(traits[5][0]), float(traits[5][1][:-1])
        closest_MJD_index = (np.abs(epochs - new_GLEP)).argmin()
        closest_MJD = epochs[closest_MJD_index]
        all_epochs = np.append(all_epochs, closest_MJD)
        size = np.abs(closest_MJD - master_traits[4])
        #print(closest_MJD)
        
        if (epoch_finding_mode == False):    
            # run tempo2 again with 0 phase MJD
            editting_par(par, closest_MJD)
            # TEMPORARY LINE - RESTRICT TO EXACT EPOCH
            #editting_par(par, 60000)
            traits = run_fit(par, tim, no_phase_fit=True, recovery_mode = True)
            print(traits)
            # traits takes the form of f0, f0_e, f1, f1_e, ph, epochs, epoch_e
            
            results = sim_arg, traits[0], traits[1], traits[2], traits[3], traits[4], num_toas, size, closest_MJD
            all_results = np.vstack((all_results, results))
            editting_par(par, 0, "GLF0_1")
            editting_par(par, 0, "GLF1_1")
            editting_par(par, 0)
            
        print(str(number+1) + ".", end="")
        sys.stdout.flush()
        
    end_time = time.time()
    print("]")
    print("done! took " + f"{(end_time - start_time):.3f} seconds")
    if (epoch_finding_mode == True): return all_epochs
    
    all_results = all_results.astype('float64')
    
    
    
    # Used this following line for quadrature error calculation but it didnt work
    # np.sqrt((np.std(all_results[:,1])/avg_f0)**2 + (np.mean(all_results[:,2]/avg_f0))**2),
                            
    return all_results

def results_averager(results):
    avg_f0 = np.mean(results[:,1])
    avg_f1 = np.mean(results[:,3])  
    # convoluted code returning average results and their errors
    avg_results = np.array([avg_f0,
                            np.std(results[:,1]),
                            avg_f1,
                            np.std(results[:,3])])
    
    return avg_results
    
    
def find_const(toas, sequence_type, const_args, sim_args, desired_toas, leeway):
    # A quick algorithm to find a constant with a given number of toas
    # We need passed args to take the form: cadence_start, offset, maxgap, const
    start_cadence, start_offset, max_gap = const_args
    min_const, max_const, toa_iterations = sim_args
    
    
    constants = np.linspace(min_const, max_const, toa_iterations)
    
    # First we must find the cadence strategy which gives a set number of TOAs
    choesn_const = 0
    given_toas = 0
    for constant in constants:
        num_toas = tim_sampling.sample_from_toas(toas, sequence_type, (start_cadence, start_offset, max_gap, constant), verbose=False, counting_mode=True)[1]
        #print(constant, num_toas)
        if num_toas < desired_toas + leeway and num_toas > desired_toas - leeway:
            choesn_const = constant
            given_toas = num_toas
            return choesn_const, given_toas
    return 0, 0

def constant_finder():
    # Code which plots out the average time between observations for a given constant, for all three of the cadence strategies  (at 20days max gap)   
    desired_abdo = 7.5
    fig = plt.figure(figsize=(16, 4))
    gs = fig.add_gridspec(1, 4, wspace=0)
    axs = gs.subplots(sharey=True)

    fig.suptitle("average days between observations for a given constant and strategy")
    fig.supylabel("average days between observations", y=0.5, x=0.1)

    # Logarithmic
    adbos = np.empty((0,1))
    constants = np.linspace(0.5, 4, 1000)
    for constant in constants:
        args = (0.5, 0, 20, constant)
        adbos = np.append(adbos, tim_sampling.fadbo('logarithmic', args))
        
    pos = np.where(np.abs(np.diff(adbos)) >= 0.5)[0]+1
    x = np.insert(constants, pos, np.nan)
    y = np.insert(adbos, pos, np.nan)    
    axs[0].plot(x, y)
    axs[0].set_xlabel("logarithmic constant")
    axs[0].set_title("logarithmic")
    axs[0].set_xlim(0.4, 4.1)
    axs[0].set_ylim(0.4, 14)
    print("log consts where adbo is 5")
    item = np.where(np.abs(y - desired_abdo) < 0.01,)
    print(x[item])
        
    # Arithmetic
    adbos = np.empty((0,1))
    constants = np.linspace(0.5, 4, 1000)
    for constant in constants:
        args = (0.5, 0, 20, constant)
        adbos = np.append(adbos, tim_sampling.fadbo('arithmetic', args))
        
    pos = np.where(np.abs(np.diff(adbos)) >= 0.5)[0]+1
    x = np.insert(constants, pos, np.nan)
    y = np.insert(adbos, pos, np.nan)    
    axs[1].plot(x, y)
    axs[1].set_xlabel("sequential increase")
    axs[1].set_title("arithmetic")
    axs[1].set_xlim(0.4, 4.1)
    print("arith consts where adbo is 5")
    item = np.where(np.abs(y - desired_abdo) < 0.01,)
    print(x[item])
    
    # geometric
    adbos = np.empty((0,1))
    constants = np.linspace(1.01, 6, 1000)
    for constant in constants:
        args = (0.5, 0, 20, constant)
        adbos = np.append(adbos, tim_sampling.fadbo('geometric', args))
        
    pos = np.where(np.abs(np.diff(adbos)) >= 0.5)[0]+1
    x = np.insert(constants, pos, np.nan)
    y = np.insert(adbos, pos, np.nan)    
    axs[2].plot(x, y)
    axs[2].set_xlabel("multiplicative increase")
    axs[2].set_title("geometric")
    axs[2].set_xlim(0.4, 6.1)
    print("geo consts where adbo is 5")
    item = np.where(np.abs(y - desired_abdo) < 0.01,)
    print(x[item])
    
    # periodic
    adbos = np.empty((0,1))
    constants = np.linspace(0.5, 20, 1000)
    for constant in constants:
        args = (0.5, 0, 20, constant)
        adbos = np.append(adbos, tim_sampling.fadbo('periodic', args))
        
    pos = np.where(np.abs(np.diff(adbos)) >= 0.5)[0]+1
    x = np.insert(constants, pos, np.nan)
    y = np.insert(adbos, pos, np.nan)    
    axs[3].plot(x, y)
    axs[3].set_xlabel("period (days)")
    axs[3].set_title("periodic")
    axs[3].set_xlim(0.4, 20.1)
    print("peri consts where adbo is 5")
    item = np.where(np.abs(y - desired_abdo) < 0.01,)
    print(x[item])
    
def diff_plot():
    # Plots our DDnu and DDnudot results for each of the cadence strategies
    
    toas = np.genfromtxt("master_toas_exp.tim", skip_header=1, usecols=[2])
    # Using pandas to read in the master file, probably a better way to do this but it works for now.
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv("master_file_exp.par", sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']))

    print(master_traits)
    
    iters = 25
    args = (0.5, 0, 20)
    
    seq = 'logarithmic'
    const = 1.0991
    passed_args = args[0], args[1], args[2], const
    print("numtoas of log", tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results = single_simulate(toas, seq, args, const, num_sps=iters)
    results = results_averager(all_results)
    plt.scatter(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], facecolors='none', edgecolors='black', s=all_results[:,7]*25, zorder=10)
    plt.errorbar(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], xerr=all_results[:,2], yerr=all_results[:,4], fmt='x', label=seq, zorder=1)    
    
    seq = 'geometric'
    const = 1.6394
    passed_args = args[0], args[1], args[2], const
    print("numtoas of "+seq, tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results = single_simulate(toas, seq, args, const, num_sps=iters)
    results = results_averager(all_results)
    plt.scatter(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], facecolors='none', edgecolors='black', s=all_results[:,7]*25, zorder=10)
    plt.errorbar(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], xerr=all_results[:,2], yerr=all_results[:,4], fmt='x', label=seq, zorder=1)    
    
    seq = 'periodic'
    const = 5
    passed_args = args[0], args[1], args[2], const
    print("numtoas of "+seq, tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results = single_simulate(toas, seq, args, const, num_sps=iters)
    results = results_averager(all_results)
    plt.scatter(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], facecolors='none', edgecolors='black', s=all_results[:,7]*25, zorder=10)
    plt.errorbar(all_results[:,1]-master_traits[0], all_results[:,3]-master_traits[1], xerr=all_results[:,2], yerr=all_results[:,4], fmt='x', label=seq, zorder=1)    
    
    plt.scatter(0, 0, c='r', label="real parameters", zorder =100)
    
    plt.xlabel(r'distance from true $\Delta \nu$')
    plt.ylabel(r'distance from true $\Delta \dot \nu$')
    plt.title(r'difference in retrieved $\Delta \nu$ and $\Delta \dot \nu$ and actual values', x=0.5, y=1.05)
    plt.legend()
    plt.savefig("figures/avg3d_nophasefit.png", dpi=400, bbox_inches="tight")
    
    
    #fig.savefig("figures/fadbos.png", dpi=400, bbox_inches="tight")
    
def histogram_plot():
    # Histogram plotter for the retrieved epochs
    toas = np.genfromtxt("master_toas.tim", skip_header=1, usecols=[2])
    numiters = 50

    fig = plt.figure(figsize=(6, 10))
    gs = fig.add_gridspec(3, 1, hspace=0)
    axs = gs.subplots(sharex = True)

    # Logarithmic
    #args = (0.5, 0, 20, 1.0991)
    results = single_simulate(toas, 'logarithmic', (0.5, 0, 20), 1.0991, num_sps=numiters, epoch_finding_mode=True)
    print(results)
    axs[0].hist(results, bins=15, range=(59999, 60006))
    axs[0].set_xlabel("epoch (MJD)")
    axs[0].set_ylabel("frequency")
    axs[0].set_title("logarithmic (const = 1.0991)")
    
    results = single_simulate(toas, 'geometric', (0.5, 0, 20), 1.6394, num_sps=numiters, epoch_finding_mode=True)
    print(results)
    axs[1].hist(results, bins=15, range=(59999, 60006))
    axs[1].set_ylabel("frequency")
    axs[1].set_title("geometric (const = 1.6394)")
    
    results = single_simulate(toas, 'periodic', (0.5, 0, 20), 5.000, num_sps=numiters, epoch_finding_mode=True)
    print(results)
    axs[2].hist(results, bins=15, range=(59999, 60006))
    axs[2].set_ylabel("frequency")
    axs[2].set_title("periodic (period = 5.000)")
    
    fig.suptitle("distributions of retrieved epochs for each cadence strategy (600 toas, 5d adbo)")
    
    fig.savefig("figures/epoch_hist_three_sharex.png", dpi=400, bbox_inches="tight")

    
def main():
    diff_plot()
    
    return

if __name__ == "__main__":
    """
    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(1, 4, wspace=0)
    axs = gs.subplots(sharey=True)
    curr = None
    """
    main()
        



