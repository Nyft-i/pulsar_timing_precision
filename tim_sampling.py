import numpy as np
import matplotlib.pyplot as plt
import time

def find_sequence_period_info(sequence_type, args):
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
        return period, 1
    else:
        print("invalid sequence type. break.")
        return 0
    
    total_time = 0
    num_toas = 0
    cadence = cadence_start
    while cadence <= max_gap:
        total_time += cadence
        num_toas += 1
        if sequence_type=='logarithmic': cadence = (np.log(1/10 * cadence + 1) * log_const) 
        elif sequence_type=='arithmetic': cadence = cadence + sequential_increase
        elif sequence_type=='exponential': cadence = np.power(cadence,exp_increase)
        elif sequence_type=='geometric': cadence = cadence * multiplicative_increase
        elif sequence_type=='periodic': cadence = period
        
    return total_time, num_toas

def fadbo(sequence_type, args):
    strat_p, num_toas = find_sequence_period_info(sequence_type, args)
    return strat_p/num_toas
    
    
def sample_from_toas(toas, sequence_type, args, verbose=False, counting_mode = False, strat_period = 0):
    # Setup
    edit_toas = toas.copy()
    end = np.max(toas) - strat_period # marker is in mjd, we want the endpoint to be one strategy period before the end of the data
    new_toas = np.zeros(0)
    num_toas = 0
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
    
    marker = np.min(edit_toas) + marker_offset 
    end = end + marker_offset # adjusts the endpoint to ensure that on average the data set is the same size
    cadence = cadence_start
    if verbose == True: print("starting cadence: " + str(cadence_start))
    #time.sleep(1)
    
    while(marker < end):
        closest_index = (np.abs(edit_toas - marker)).argmin()
        
        # Checks if the closest index has already been picked before.
        if(counting_mode==False):
            if((np.isin(closest_index, indexes))):
                # Currently just skips over it if so, could implement to find the next closest.
                if verbose==True: print("double counted! skipping")
            else:
                # Appends that particular TOA to the new list of empty ToAs.
                new_toas = np.append(new_toas, edit_toas[closest_index])
                # Ads the index also to avoid double counting
                indexes = np.append(indexes, closest_index)
                # Removes the ToA from the list so it cant be picked again, does this by setting its value to infinity so it is never picked again.
                edit_toas[closest_index] = float("inf")
                cadence_list = np.append(cadence_list, cadence)
        
        if sequence_type=='logarithmic': cadence = (np.log(1/10 * cadence + 1) * log_const) 
        elif sequence_type=='arithmetic': cadence = cadence + sequential_increase
        elif sequence_type=='geometric': cadence = cadence * multiplicative_increase
        elif sequence_type=='periodic': cadence = period
        
        if(cadence > max_gap): cadence = cadence_start
        if verbose==True: print("current cadence: " + str(cadence))
        marker += cadence
        num_toas += 1
        #if verbose: time.sleep(0.5)
                
    return indexes, num_toas

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
    return len(indexes)

#print("i do this!")
def main():
    # takes user input for sampling
    timfile = input("Enter the name of the tim file you wish to sample: ")
    sequence_type = input("Enter the type of sequence you wish to sample with: ")
    cadence_start = float(input("Enter the starting cadence: "))
    marker_offset = float(input("Enter the marker offset: "))
    max_gap = float(input("Enter the maximum gap: "))
    if sequence_type == 'logarithmic':
        log_const = float(input("Enter the logarithmic constant: "))
        args = [cadence_start, marker_offset, max_gap, log_const]
    elif sequence_type == 'arithmetic':
        sequential_increase = float(input("Enter the sequential increase: "))
        args = [cadence_start, marker_offset, max_gap, sequential_increase]
    elif sequence_type == 'geometric':
        multiplicative_increase = float(input("Enter the multiplicative increase: "))
        args = [cadence_start, marker_offset, max_gap, multiplicative_increase]
    elif sequence_type == 'periodic':
        period = float(input("Enter the period: "))
        args = [cadence_start, marker_offset, max_gap, period]
    else:
        print("invalid sequence type. break.")
        return 0
    
    # Reads the .tim file and extracts the TOAs
    toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])
    indexes, num_toas = sample_from_toas(toas, sequence_type, args, verbose=True)
    gen_new_tim(timfile, indexes, "new.tim")
    

if __name__ == "__main__":
    main()