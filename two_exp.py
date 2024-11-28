# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:08:57 2024

@author: laure
"""

import tim_sampling
import numpy as np
import subprocess
import pandas
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import os
import time

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
    #midpoint code 
    
    toas = np.genfromtxt(tim, skip_header=1, usecols=[2])
    difference = toas - master_traits[3]
    
    first_post_i = np.abs(np.where(difference>0, difference, np.inf)).argmin()
    last_pre_i = np.abs(np.where(difference<0, difference, -np.inf)).argmin()
    
    first_post = toas[first_post_i]
    last_pre = toas[last_pre_i]
    
    mid_point = (first_post + last_pre)/2
    
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
    
def run_fit(par, tim, recovery_mode = False):
    command = [
        "tempo2",
        "-f", par, tim,
        "-nofit",
        "-fit", "GLF0_1",
        "-fit", "GLF1_1",
        "-fit", "GLPH_1",
        "-fit", "F1",
        "-fit", "F0"]
    
    if recovery_mode == True :
        command_rec = ["-fit", "GLF0D_1",
                       "-fit", "GLTD_1",
                       "-fit", "GLF0D_2",
                       "-fit", "GLTD_2"]
        command = np.hstack((command,command_rec))
            

    #print(' '.join(command), file=sys.stderr)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    out, err = proc.communicate()
    all_fields = out.split("\n")    
    #print(command)

    f0, f0_e, f1, f1_e, ph, epochs, epoch_e, recovered_F0, recovered_F0_e, recovered_timescale, recovered_timescale_e, pulsar_f0, pulsar_f1, recovered_F0_2, recovered_F0_2_e, recovered_timescale_2, recovered_timescale_2_e = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0                        

    for this_field in all_fields:
        fields = this_field.split()
        
        if len(fields) > 2.0:
            if fields[0] == "GLF0_1" and f0 == 0:
                f0 = fields[2]
                f0_e = fields[3]
            if fields[0] == "GLF1_1" and f1 == 0:
                f1 = fields[2]
                f1_e = fields[3]
            if fields[0] == "GLPH_1" and ph == 0:
                ph = fields[2]
            if fields[0] == "MJD":
                epochs = fields[7], fields[9]
                epoch_e = fields[12]
            if fields[0] == "GLF0D_1" and recovered_F0 == 0:
                recovered_F0 = fields[2]
                recovered_F0_e = fields[3]
            if fields[0] == "GLTD_1" and recovered_timescale == 0:
                recovered_timescale = fields[2]
                recovered_timescale_e = fields[3]
            if fields[0] == "F0" and pulsar_f0 == 0:
                pulsar_f0 = fields[3]
                pulsar_f0_e = fields[4]
            if fields[0] == "F1" and pulsar_f1 == 0:
                pulsar_f1 = fields[3]
                pulsar_f1_e = fields[4]
            if fields[0] == "GLF0D_2" and recovered_F0_2 == 0:
                recovered_F0_2 = fields[2]
                recovered_F0_2_e = fields[3]
            if fields[0] == "GLTD_2" and recovered_timescale_2 == 0:
                recovered_timescale_2 = fields[2]
                recovered_timescale_2_e = fields[3]        
        
    #print(f0, f0_e, f1, f1_e, ph, epochs, epoch_e, recovered_F0, recovered_F0_e, recovered_timescale, recovered_timescale_e, pulsar_f0, pulsar_f1, recovered_F0_2, recovered_F0_2_e, recovered_timescale_2, recovered_timescale_2_e)
    try:
        if recovery_mode == True:
            return f0, f0_e, f1, f1_e, ph, epochs, epoch_e, recovered_F0, recovered_F0_e, recovered_timescale, recovered_timescale_e, pulsar_f0, pulsar_f0_e, pulsar_f1, pulsar_f1_e, recovered_F0_2, recovered_F0_2_e, recovered_timescale_2, recovered_timescale_2_e

        return f0, f0_e, f1, f1_e, ph, epochs, epoch_e,0,0,0,0, pulsar_f0, pulsar_f0_e, pulsar_f1, pulsar_f1_e,0,0,0,0
    except UnboundLocalError:
        return None

def single_simulate(toas, sequence_type, const_args, sim_arg, recovery, verbose = False, master_tim="master_toas_2exp.tim", master_par="glitchC_master.par", temp_par = "glitchC_temp.par", num_sps=1, epoch_finding_mode=False):
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
    #print(master_traits)
    # adds some 5d random variation so that we dont run into issues with the sample being the same every time
    passed_args = const_args[0], const_args[1], const_args[2], sim_arg
    strategy_period, strat_toas = tim_sampling.find_sequence_period_info(sequence_type, passed_args)
    start_randomiser = np.random.randint(0, strategy_period*100, (num_sps))
    start_randomiser = start_randomiser/100
    all_results = np.zeros((0,21))
    all_epochs = np.zeros(0)
    
    print("[",end="")
    # For each offset, we generate a new set of toas, run tempo2, and compare the results to the master file
    #print("running simulation for "+sequence_type+" sequence type\n[",end="")
    for number, offset in enumerate(start_randomiser):
        print(str(number+1)+".",end="")
        # force python console to update
        sys.stdout.flush()
        # We need passed args to take the form: cadence_start, offset, maxgap, const
        # const_args: start cadence, start offset, max_gap
        passed_args = const_args[0], const_args[1]+offset, const_args[2], sim_arg
        
        indexes, num_toas = tim_sampling.sample_from_toas(toas, sequence_type, passed_args, verbose, strat_period=strategy_period)
        
        temp_tim = sequence_type+"_toas.tim"
        #print(indexes)
        tim_sampling.gen_new_tim(master_tim, indexes, temp_tim)
        
        par, tim = temp_par, temp_tim
        # Ensure the par file is clean
        editting_par(par, 0, "GLF0_1")
        editting_par(par, 0, "GLF1_1")
        editting_par(par, 0)
        editting_par(par, 0, "GLEP_2")
        
        # Residual loading glep finder code, put it in the par file
        new_GLEP = epoch_finder(par, tim, master_traits)
        #all_epochs = np.append(all_epochs, new_GLEP)        
        #print(new_GLEP)
        editting_par(par, new_GLEP)
        editting_par(par, new_GLEP, "GLEP_2")
        
        # run tempo2
        traits = run_fit(par, tim)
        #print(traits)
        editting_par(par, traits[0], "GLF0_1")
        editting_par(par, traits[2], "GLF1_1")
        
        
        epochs = float(traits[5][0]), float(traits[5][1][:-1])
        closest_MJD_index = (np.abs(epochs - new_GLEP)).argmin()
        closest_MJD = epochs[closest_MJD_index]
        all_epochs = np.append(all_epochs, closest_MJD)
        size = np.abs(closest_MJD - master_traits[4])
        #print("mjd used:",closest_MJD)
        
        if (epoch_finding_mode == False):    
            # run tempo2 again with 0 phase MJD
            editting_par(par, closest_MJD)
            editting_par(par, closest_MJD, "GLEP_2")
            traits = run_fit(par, tim, recovery_mode = recovery)
            #print(traits)
            # traits takes the form of f0, f0_e, f1, f1_e, ph, epochs, epoch_e
            # results takes the form sim_arg, df0, df0e, df1, df1e, phase, numtoas, size, closestmjd, recoveryf0, recoveryf0e, recoveryt, recoveryte
            results = sim_arg, traits[0], traits[1], traits[2], traits[3], traits[4], num_toas, size, closest_MJD, traits[7], traits[8], traits[9], traits[10], traits[11], traits[12], traits[13], traits[14], traits[15], traits[16], traits[17], traits[18]
            all_results = np.vstack((all_results, results))
        
        # clean up at the end also
        editting_par(par, 0, "GLF0_1")
        editting_par(par, 0, "GLF1_1")
        editting_par(par, 0)
        editting_par(par, 0, "GLEP_2")
            
        #print("finished pulsar ", number)
    print("]")
    end_time = time.time()
    #print("]")
    print("done! took " + f"{(end_time - start_time):.3f} seconds")
    if (epoch_finding_mode == True): return all_epochs
    
    all_results = all_results.astype('float64')
    
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

def diff_plot_recoveries():
    # Plots our DDnu and DDnudot results for each of the cadence strategies
    
    par = "glitchC_master.par"
    tim = "master_toas_2exp.tim"
    temppar = "glitchC_temp.par"
    toas = np.genfromtxt(tim, skip_header=1, usecols=[2])
    # Using pandas to read in the master file, probably a better way to do this but it works for now.
    cols = ["Element Name", "Value", "Fitting", "Error"]
    master_properties = pandas.read_csv(par, sep="\s+", names=cols)
    master_traits = (float(master_properties.loc[master_properties['Element Name'] == "GLF0_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLF1_1"]['Value']), 
                    float(master_properties.loc[master_properties['Element Name'] == "GLPH_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "PEPOCH"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLEP_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_1"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLF0D_2"]['Value']),
                    float(master_properties.loc[master_properties['Element Name'] == "GLTD_2"]['Value']))
    print(master_traits)
    
    iters = 50
    args = (0.5, 0, 20)
    
    fig = plt.figure(figsize=(9, 3))
    gs = fig.add_gridspec(1, 3, wspace = 0)
    axs = gs.subplots(sharey = True, sharex = True)
    
    fig.suptitle(r'difference in retrieved recovery portion of $\Delta \nu$ and $\tau_r$ and actual values', x=0.5, y=1.05)
    fig.supylabel(r'$\Delta \nu_d$ - ' + str(master_traits[5]), y=0.45, x=0.06)
    fig.supxlabel(r'$\tau_r$ - ' + str(master_traits[6]), y = -0.05)
    
    seq = 'logarithmic'
    const = 25.7197
    passed_args = args[0], args[1], args[2], const
    print("numtoas of log", tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results_log = single_simulate(toas, seq, args, const, True, num_sps=iters, master_par=par, master_tim=tim, temp_par = temppar)
    x_first_exp_log = np.mean(all_results_log[:,11]) - master_traits[6]
    y_first_exp_log = np.mean(all_results_log[:,9]) - master_traits[5]
    
    x_second_exp_log = np.mean(all_results_log[:,19]) - master_traits[8]
    y_second_exp_log = np.mean(all_results_log[:,17]) - master_traits[7]
    
    x_err_first_log = np.std(all_results_log[:,11])
    y_err_first_log = np.std(all_results_log[:,9])
    
    x_err_second_log = np.std(all_results_log[:,19])
    y_err_second_log = np.std(all_results_log[:,17])
    
    results_log = results_averager(all_results_log)
    
    
    
    # df0 and df1
    axs[0].scatter(all_results_log[:,11]-master_traits[6], all_results_log[:,9]-master_traits[5], facecolors='none', edgecolors='mediumorchid', s=all_results_log[:,7]*25, zorder=10, alpha = 0.3)
    axs[0].errorbar(all_results_log[:,11]-master_traits[6], all_results_log[:,9]-master_traits[5], xerr=all_results_log[:,12], yerr=all_results_log[:,10], fmt='x', zorder=1, alpha = 0.3, color = "mediumorchid")    
    axs[0].errorbar(x_first_exp_log, y_first_exp_log, xerr = x_err_first_log, yerr = y_err_first_log, zorder = 50, fmt = "x", color = "darkmagenta", label = "mean")

    axs[0].set_title("logarithmic")
    
    #seq = 'arithmetic'
    #const = 2.3744
    #passed_args = args[0], args[1], args[2], const
    #print("numtoas of log", tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    #all_results_arith = single_simulate(toas, seq, args, const, True, num_sps=iters, master_par=par, master_tim=tim, temp_par = temppar)
    #x_avg = np.mean(all_results_arith[:,11]) - master_traits[6]
    #y_avg = np.mean(all_results_arith[:,9]) - master_traits[5]
    #x_err = np.std(all_results_arith[:,11])
    #y_err = np.std(all_results_arith[:,9])
    
    #results_arith = results_averager(all_results_arith)
    
    # df0 and df1
    #axs[0].scatter(all_results_arith[:,11]-master_traits[6], all_results_arith[:,9]-master_traits[5], facecolors='none', edgecolors='tab:blue', s=all_results_arith[:,7]*25, zorder=10, alpha = 0.3)
    #axs[0].errorbar(all_results_arith[:,11]-master_traits[6], all_results_arith[:,9]-master_traits[5], xerr=all_results_arith[:,12], yerr=all_results_arith[:,10], fmt='x', label=seq, zorder=1, alpha = 0.3, color = "tab:blue")    
    #axs[0].errorbar(x_avg, y_avg, xerr = x_err, yerr = y_err, label = seq, zorder = 50, fmt = "x", color = "darkblue")
    
    #axs[0].set_title("arithmetic")
    
    # timescale on x and recovery 
    seq = 'geometric'
    const = 1.6394
    passed_args = args[0], args[1], args[2], const
    print("numtoas of "+seq, tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results_geo = single_simulate(toas, seq, args, const, True, num_sps=iters, master_par=par, master_tim=tim, temp_par = temppar)
    x_first_exp_geo = np.mean(all_results_geo[:,11]) - master_traits[6]
    y_first_exp_geo = np.mean(all_results_geo[:,9]) - master_traits[5]
    
    print(all_results_geo[:,11])
    print(np.mean(all_results_geo[:,11]))
    
    print(all_results_geo[:,9])
    print(np.mean(all_results_geo[:,9]))
    
    
    x_second_exp_geo = np.mean(all_results_geo[:,19]) - master_traits[8]
    y_second_exp_geo = np.mean(all_results_geo[:,17]) - master_traits[7]
    
    x_err_second_geo = np.std(all_results_geo[:,19])
    y_err_second_geo = np.std(all_results_geo[:,17])
    
    x_err_first_geo = np.std(all_results_geo[:,11])
    y_err_first_geo = np.std(all_results_geo[:,9])
    
    results_geo = results_averager(all_results_geo)
   
    # df0 and df1
    axs[2].scatter(all_results_geo[:,11]-master_traits[6], all_results_geo[:,9]-master_traits[5],  facecolors='none', edgecolors='orange', s=all_results_geo[:,7]*25, zorder=10, alpha = 0.3)
    axs[2].errorbar(all_results_geo[:,11]-master_traits[6], all_results_geo[:,9]-master_traits[5], xerr=all_results_geo[:,12], yerr=all_results_geo[:,10], fmt='x', zorder=1, alpha = 0.3, color = "orange")    
    axs[2].errorbar(x_first_exp_geo, y_first_exp_geo, xerr = x_err_first_geo, yerr = y_err_first_geo, zorder = 50, fmt = "x", color = "goldenrod")

    axs[2].set_title("geometric")
    
    seq = 'periodic'
    const = 5
    passed_args = args[0], args[1], args[2], const
    print("numtoas of "+seq, tim_sampling.sample_from_toas(toas, seq, passed_args, counting_mode=True)[1])
    all_results_per = single_simulate(toas, seq, args, const, True, num_sps=iters, master_par=par, master_tim=tim, temp_par = temppar)
    x_first_exp_per = np.mean(all_results_per[:,11]) - master_traits[6]
    y_first_exp_per = np.mean(all_results_per[:,9]) - master_traits[5]
    
    x_second_exp_per = np.mean(all_results_per[:,19]) - master_traits[8]
    y_second_exp_per = np.mean(all_results_per[:,17]) - master_traits[7]
    
    x_err_second_per = np.std(all_results_per[:,19])
    y_err_second_per = np.std(all_results_per[:,17])
    
    x_err_first_per = np.std(all_results_per[:,11])
    y_err_first_per = np.std(all_results_per[:,9])
    
    results_per = results_averager(all_results_per)
    
    # df0 and df1
    axs[1].scatter(all_results_per[:,11]-master_traits[6], all_results_per[:,9]-master_traits[5],  facecolors='none', edgecolors='limegreen', s=all_results_per[:,7]*25, zorder=10, alpha = 0.3)
    axs[1].errorbar(all_results_per[:,11]-master_traits[6], all_results_per[:,9]-master_traits[5], xerr=all_results_per[:,12], yerr=all_results_per[:,10], fmt='x', zorder=1, alpha = 0.3, color = "limegreen")    
    axs[1].errorbar(x_first_exp_per, y_first_exp_per, xerr = x_err_first_per, yerr = y_err_first_per, zorder = 50, fmt = "x", color = "darkgreen")
    
    axs[0].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[1].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[2].scatter(0, 0, c='r', label="real parameters", zorder =100)
    
    axs[1].set_title("periodic")
    
    #axs[2].legend()
    axs[0].legend()
    
    plt.savefig("figures/2recovery_params_3d_w_average.png", dpi=400, bbox_inches="tight") 
    
    plt.clf()
    
    fig = plt.figure(figsize=(9, 3))
    gs = fig.add_gridspec(1, 3, wspace = 0)
    axs = gs.subplots(sharey = True, sharex = True)
    
    fig.suptitle(r'difference in retrieved $\Delta \nu$ and $\Delta \dot \nu$ and actual values', x=0.5, y=1.05)
    fig.supylabel(r'$\Delta \dot \nu - $' + str(master_traits[1]), y=0.5, x=0.06)
    fig.supxlabel(r'$\Delta \nu - $' + str(master_traits[0]), y = -0.13)
    
    axs[0].scatter(all_results_log[:,1]-master_traits[0], all_results_log[:,3]-master_traits[1], facecolors='none', edgecolors='mediumorchid', s=all_results_log[:,7]*25, zorder=10, alpha = 0.3)
    axs[0].errorbar(all_results_log[:,1]-master_traits[0], all_results_log[:,3]-master_traits[1], xerr=all_results_log[:,2], yerr=all_results_log[:,4], fmt='x', zorder=1, alpha = 0.3, color = "mediumorchid")    
    axs[0].errorbar(results_log[0]-master_traits[0], results_log[2] - master_traits[1], xerr = results_log[1], yerr = results_log[3], zorder = 50, fmt = "x", color = "darkmagenta", label = "mean")
   
    #axs[0].scatter(all_results_arith[:,1]-master_traits[0], all_results_arith[:,3]-master_traits[1], facecolors='none', edgecolors='tab:blue', s=all_results_arith[:,7]*25, zorder=10, alpha = 0.3)
    #axs[0].errorbar(all_results_arith[:,1]-master_traits[0], all_results_arith[:,3]-master_traits[1], xerr=all_results_arith[:,2], yerr=all_results_arith[:,4], fmt='x', label=seq, zorder=1, alpha = 0.3, color = "tab:blue")    
    #axs[0].errorbar(results_arith[0]-master_traits[0], results_arith[2] - master_traits[1], xerr = results_arith[1], yerr = results_arith[3],label = seq, zorder = 50, fmt = "x", color = "darkblue")
    
    axs[2].scatter(all_results_geo[:,1]-master_traits[0], all_results_geo[:,3]-master_traits[1], facecolors='none', edgecolors='orange', s=all_results_geo[:,7]*25, zorder=10, alpha = 0.3)
    axs[2].errorbar(all_results_geo[:,1]-master_traits[0], all_results_geo[:,3]-master_traits[1], xerr=all_results_geo[:,2], yerr=all_results_geo[:,4], fmt='x', zorder=1, alpha = 0.3, color = "orange")
    axs[2].errorbar(results_geo[0]-master_traits[0], results_geo[2] - master_traits[1], xerr = results_geo[1], yerr = results_geo[3],zorder = 50, fmt = "x", color = "goldenrod")
   
    axs[1].scatter(all_results_per[:,1]-master_traits[0], all_results_per[:,3]-master_traits[1], facecolors='none', edgecolors='limegreen', s=all_results_per[:,7]*25, zorder=10, alpha = 0.3)
    axs[1].errorbar(all_results_per[:,1]-master_traits[0], all_results_per[:,3]-master_traits[1], xerr=all_results_per[:,2], yerr=all_results_per[:,4], fmt='x', zorder=1, alpha = 0.3, color = "limegreen")
    axs[1].errorbar(results_per[0]-master_traits[0], results_per[2] - master_traits[1], xerr = results_per[1], yerr = results_per[3], zorder = 50, fmt = "x", color = "darkgreen")
 
    axs[0].set_title("logarithmic")
    #axs[0].set_title("arithmetic")
    axs[2].set_title("geometric")
    axs[1].set_title("periodic")
    
    axs[0].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[1].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[2].scatter(0, 0, c='r', label="real parameters", zorder =100)
    
    #axs[2].legend()
    axs[0].legend()
    
    plt.savefig("figures/2recovery_normal_params_3d_w_average.png", dpi=400, bbox_inches="tight")
    
    plt.clf()
    
    fig = plt.figure(figsize=(9, 3))
    gs = fig.add_gridspec(1, 3, wspace = 0)
    axs = gs.subplots(sharey = True, sharex = True)
    
    fig.suptitle(r'difference in retrieved recovery portion of $\Delta \nu$ and $\tau_r$ and actual values for long recovery portion', x=0.5, y=1.05)
    fig.supylabel(r'$\Delta \nu_d$ - ' + str(master_traits[7]), y=0.45, x=0.06)
    fig.supxlabel(r'$\tau_r$ - ' + str(master_traits[8]), y = -0.05)
    
    axs[0].scatter(all_results_log[:,19]-master_traits[8], all_results_log[:,17]-master_traits[7], facecolors='none', edgecolors='mediumorchid', s=all_results_log[:,7]*25, zorder=10, alpha = 0.3)
    axs[0].errorbar(all_results_log[:,19]-master_traits[8], all_results_log[:,17]-master_traits[7], xerr=all_results_log[:,20], yerr=all_results_log[:,18], fmt='x', zorder=1, alpha = 0.3, color = "mediumorchid")    
    axs[0].errorbar(x_second_exp_log, y_second_exp_log, xerr = x_err_second_log, yerr = y_err_second_log, zorder = 50, fmt = "x", color = "darkmagenta", label = "mean")

    axs[0].set_title("logarithmic")
    
    axs[2].scatter(all_results_geo[:,19]-master_traits[8], all_results_geo[:,17]-master_traits[7],  facecolors='none', edgecolors='orange', s=all_results_geo[:,7]*25, zorder=10, alpha = 0.3)
    axs[2].errorbar(all_results_geo[:,19]-master_traits[8], all_results_geo[:,17]-master_traits[7], xerr=all_results_geo[:,20], yerr=all_results_geo[:,18], fmt='x', zorder=1, alpha = 0.3, color = "orange")    
    axs[2].errorbar(x_second_exp_geo, y_second_exp_geo, xerr = x_err_second_geo, yerr = y_err_second_geo, zorder = 50, fmt = "x", color = "goldenrod")

    axs[2].set_title("geometric")
    
    axs[1].set_title("periodic")
    
    axs[1].scatter(all_results_per[:,19]-master_traits[8], all_results_per[:,17]-master_traits[7],  facecolors='none', edgecolors='limegreen', s=all_results_per[:,7]*25, zorder=10, alpha = 0.3)
    axs[1].errorbar(all_results_per[:,19]-master_traits[8], all_results_per[:,17]-master_traits[7], xerr=all_results_per[:,20], yerr=all_results_per[:,18], fmt='x', zorder=1, alpha = 0.3, color = "limegreen")    
    axs[1].errorbar(x_second_exp_per, y_second_exp_per, xerr = x_err_second_per, yerr = y_err_second_per, zorder = 50, fmt = "x", color = "darkgreen")
    
    axs[0].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[1].scatter(0, 0, c='r', label="real parameters", zorder =100)
    axs[2].scatter(0, 0, c='r', label="real parameters", zorder =100)
    
    #axs[2].legend()
    axs[0].legend()
    
    plt.savefig("figures/2recovery_second_params_3d_w_average.png", dpi=400, bbox_inches="tight") 

def data_output():
    #simulation params
    seq = "periodic"
    tim_iters = 230
    sub_iters = 45
    const = 15
    max_gap = 30
    start_cad = 2
    
    #glitch params
    tim_name = "master_toas_2exp.tim"
    par_file = "glitchC_master.par"
    temp_file = "glitchC_temp.par"
    
    #other params
    args = (start_cad, 0, max_gap)
    par_file_no_fileext = par_file.split(".")[0]
    total_sims = tim_iters*sub_iters
    curr_time = time.strftime("%H:%M")
    old_name = seq+"_"+str(const)+"_"+par_file_no_fileext+"_"+str(total_sims)+"s_"+str(curr_time)+".txt"

    # note start time
    start_time = time.time()
    # loop over all values lower than tim_iters
    for curr_tim in range(tim_iters):
        tim_sampling.gen_fresh_toas(par_file, tim_name)
        print("tim file '"+tim_name+"' generated")
        toas = np.genfromtxt(tim_name, skip_header=1, usecols=[2])
        print("starting sub-simulations for tim file "+str(curr_tim+1)+".")
        all_results = single_simulate(toas, seq, args, const, True, num_sps = sub_iters, temp_par=temp_file, master_par=par_file, master_tim=tim_name)
        print("finished sub-simulations for tim file "+str(curr_tim+1)+".")
        f=open(old_name,'a')
        np.savetxt(f, all_results, fmt = "%s", delimiter = " ")
        f.close()
        curr_time = time.strftime("%H:%M")
        new_name = seq+"_"+str(const)+"_"+par_file_no_fileext+"_"+str(total_sims)+"s_"+str(curr_time)+".txt"
        os.rename(old_name, new_name)
        print("appended data to file: "+new_name)
        old_name = new_name
    
    # note end time
    end_time = time.time()
    time_diff = end_time - start_time
    # represetnt time taken in hours, minutes and seconds
    hours = time_diff//3600
    minutes = (time_diff%3600)//60
    seconds = time_diff%60
    
    print("all simulations complete! total time taken: " +f"{hours:.0f}h {minutes:.0f}m {seconds:.0f}s")
          
    #print the time to one 
    print()       

def main():
    data_output()

    
    return

main()







