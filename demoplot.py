import matplotlib.pyplot as plt
import numpy as np
import time
import tim_sampling

def event_maker(sequence_type, args):
    if sequence_type == 'logarithmic':
        cadence_start, marker_offset, max_gap, log_const = args
    elif sequence_type == 'arithmetic':
        cadence_start, marker_offset, max_gap, sequential_increase = args
    elif sequence_type == 'exponential':
        cadence_start, marker_offset, max_gap, exp_increase = args
    elif sequence_type == 'geometric':
        cadence_start, marker_offset, max_gap, multiplicative_increase = args
    elif sequence_type == 'periodic':
        cadence_start, marker_offset, max_gap, period = args
    else:
        print("invalid sequence type. break.")
        
    start = marker_offset
    end = 150
    cadence = cadence_start

    all_obs = np.empty(0)
    endpoints = np.empty(0)

    marker = start 
    while marker <= end:
        all_obs = np.append(all_obs, marker)
        #print("current cadence: ",cadence)

        if sequence_type=='logarithmic': cadence = (np.log(1/10 * cadence + 1) * log_const) 
        elif sequence_type=='arithmetic': cadence = cadence + sequential_increase
        elif sequence_type=='exponential': cadence = np.power(cadence,exp_increase)
        elif sequence_type=='geometric': cadence = cadence * multiplicative_increase
        elif sequence_type=='periodic': cadence = period
        
        if(cadence > max_gap): 
            endpoints = np.append(endpoints, marker)
            cadence = cadence_start

        marker += cadence
        
    return all_obs, endpoints

fig = plt.figure(figsize=(7, 3))
gs = fig.add_gridspec(3, hspace=1)
axs = gs.subplots(sharex = True)

args = (0.5, 0, 10, 1.5)
seq = 'arithmetic'
print(tim_sampling.find_sequence_period_info(seq, args))


    
args = (0.5, 0, 20, 5)
seq = 'periodic'
print(tim_sampling.find_sequence_period_info(seq, args))

obs, ends = event_maker(seq, args)
axs[0].set_xlim(0,140)
axs[0].set_ylim(0.5,1)
axs[0].set_yticks([])
axs[0].set_title(seq, loc="left")
axs[0].eventplot(obs, linewidths=1, colors = "black", alpha=0.3)
axs[0].eventplot(ends, colors="green", zorder=5)

args = (0.5, 0, 20, 25.7197)
seq = 'logarithmic'
print(tim_sampling.find_sequence_period_info(seq, args))
obs, ends = event_maker(seq, args)
toofar = obs[12] + 21.85
axs[1].set_xlim(0,140)
axs[1].set_ylim(0.5,1)
axs[1].set_yticks([])
axs[1].set_title(seq, loc="left")
axs[1].eventplot(obs, linewidths=1, colors = "black", alpha=0.3)
axs[1].axvline(x = toofar, linewidth=1, c = "red")
axs[1].eventplot(ends, colors="green", zorder=5)

args = (0.5, 0, 20, 1.6394)
seq = 'geometric'
print(tim_sampling.find_sequence_period_info(seq, args))

obs, ends = event_maker(seq, args)
axs[2].set_xlim(0,140)
axs[2].set_ylim(0.5,1)
axs[2].set_yticks([])
axs[2].set_title(seq, loc="left")
axs[2].eventplot(obs, linewidths=1, colors = "black", zorder=1, alpha=0.3)
axs[2].eventplot(ends, colors="green", zorder=5)

fig.supxlabel("days", y=-0.05)

plt.savefig("demoplot.png", dpi=400, bbox_inches="tight")
plt.show()

        