import numpy as np
import matplotlib.pyplot as plt

def sample_from_toas(toas, sequence_type, args, verbose=False):
    # Setup
    end = np.max(toas)
    new_toas = np.zeros(0)
    indexes = np.zeros(0,dtype=int)
    cadence_list = np.zeros(0)
    
    if sequence_type == 'logarithmic':
        cadence_start, marker_offset, max_gap, log_const = args
    elif sequence_type == 'arithmetic':
        cadence_start, marker_offset, max_gap, sequential_increase = args
    elif sequence_type == 'geometric':
        cadence_start, marker_offset, max_gap, multiplicative_increase = args
    elif sequence_type == 'periodic':
        cadence_start, marker_offset, max_gap, period = args
    else:
        print("invalid sequence type. break.")
        return indexes
    
    marker = np.min(toas) + marker_offset
    cadence = cadence_start
    
    while(marker < end):
        closest_index = (np.abs(toas - marker)).argmin()
        
        # Checks if the closest index has already been picked before.
        if((np.isin(closest_index, indexes))):
            # Currently just skips over it if so, could implement to find the next closest.
            if verbose==True: print("double counted! skipping")
        else:
            # Appends that particular TOA to the new list of empty ToAs.
            new_toas = np.append(new_toas, toas[closest_index])
            # Ads the index also to avoid double counting
            indexes = np.append(indexes, closest_index)
            # Removes the ToA from the list so it cant be picked again, does this by setting its value to infinity so it is never picked again.
            toas[closest_index] = float("inf")
            cadence_list = np.append(cadence_list, cadence)
        
        if sequence_type=='logarithmic': cadence = np.exp(np.log(cadence) + log_const)
        elif sequence_type=='arithmetic': cadence = cadence + sequential_increase
        elif sequence_type=='geometric': cadence = cadence * multiplicative_increase
        elif sequence_type=='periodic': cadence = period
        
        if(cadence > max_gap): cadence = cadence_start
        if verbose==True: print("current cadence: " + str(cadence))
        marker += cadence
    return indexes

def gen_new_tim(timfile, indexes, newfile):
    # Creates a string array with identical formatting to the input .tim file.
    lines = np.genfromtxt(timfile, skip_header=1, delimiter="no-delim", dtype=str)
    new_lines = lines[indexes]

    # Re-ads the header to the top of the .tim file - unsure if this is important
    header = "FORMAT 1"
    new_lines = np.hstack((header, new_lines))

    # Puts the new tim file array into an actual tim file, ready to be read by tempo2.
    #file_name = str(start_cadence) + "_day_("+str(marker_offset)+"_offset).tim"
    np.savetxt(newfile, new_lines, fmt="%s")


timfile = "master_toas.tim"
toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])

# User changeable 
# 'logarithmic', 'arithmetic', 'geometric', 'periodic'
SEQUENCE_TYPE = 'periodic'
cadence_start = 0.5
marker_offset = 0
max_gap = 20
verbose = True

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
    indexes = sample_from_toas(toas, 'logarithmic', args, verbose)
elif SEQUENCE_TYPE == 'arithmetic':
    args = (cadence_start, marker_offset, max_gap, sequential_increase)
    indexes = sample_from_toas(toas, 'arithmetic', args, verbose)
elif SEQUENCE_TYPE == 'geometric':
    args = (cadence_start, marker_offset, max_gap, multiplicative_increase)
    indexes = sample_from_toas(toas, 'geometric', args, verbose)
elif SEQUENCE_TYPE == 'periodic':
    args = (cadence_start, marker_offset, max_gap, period)
    indexes = sample_from_toas(toas, 'periodic', args, verbose)
else:
    print("invalid sequence type. doing nothing.")    


print("number of toas: " + str(len(indexes)))

new_filename = SEQUENCE_TYPE + "_toas.tim"
gen_new_tim(timfile, indexes, new_filename)

