import tkinter as tk
from tkinter import ttk
#import threading
import numpy as np

#import simulating_cadence

app = tk.Tk()
app.title("Pulsar Timer")
app.geometry("400x500")

def forget_all():
    log_frame.grid_forget()
    arith_frame.grid_forget()
    geom_frame.grid_forget()
    period_frame.grid_forget()
    default_frame.grid_forget()

def log_frame_command():
    forget_all()
    log_frame.grid(column=1, row=0, sticky="ew")
    log_frame.grid_propagate(False)

def arith_frame_command():
    forget_all()
    arith_frame.grid(column=1, row=0, sticky="ew")
    arith_frame.grid_propagate(False)

def geom_frame_command():
    forget_all()
    geom_frame.grid(column=1, row=0, sticky="ew")
    geom_frame.grid_propagate(False)

def period_frame_command():
    forget_all()
    period_frame.grid(column=1, row=0, sticky="ew")
    period_frame.grid_propagate(False)

def simulate():
    """
    timfile = "master_toas_2.tim"
    toas = np.genfromtxt(timfile, skip_header=1, usecols=[2])



    const_args = (float(start_cadence.get()), float(start_offset.get()), float(max_gap.get()))
    sim_args = (float(log_const_min.get()), float(log_const_max.get()), int(num_iterations.get()))
    print("simulating with following args:")
    print("const_args: ")
    print(const_args)
    print("sim_args: ")
    print(sim_args)

    app.quit()
    simulating_cadence.simulate(toas, SEQUENCE_TYPE.get(), const_args, sim_args, simprogress)
    """

frame_radio = tk.Frame(app, width=200, height=100)
frame_radio.grid(column=0, row=0,sticky="ew")
frame_radio.grid_propagate(False)
frame_radio.columnconfigure(0, weight=1)

default_frame = tk.Frame(app, width=200, height=100)
default_frame.grid(column=1, row=0, sticky="ew")
default_frame.grid_propagate(False)



# Create 4 radio buttons in the left frame
SEQUENCE_TYPE = tk.StringVar()
SEQUENCE_TYPE.set('none')
ttk.Radiobutton(frame_radio, text='Logarithmic', variable=SEQUENCE_TYPE, value='logarithmic',command=log_frame_command).grid(row=0)
ttk.Radiobutton(frame_radio, text='Arithmetic', variable=SEQUENCE_TYPE, value='arithmetic',command=arith_frame_command).grid(row=1)
ttk.Radiobutton(frame_radio, text='Geometric', variable=SEQUENCE_TYPE, value='geometric',command=geom_frame_command).grid(row=2)
ttk.Radiobutton(frame_radio, text='Periodic', variable=SEQUENCE_TYPE, value='periodic',command=period_frame_command).grid(row=3)

# Log Settings
log_frame = tk.Frame(app, width=200, height=100)
log_frame.grid_propagate(False)
log_frame.columnconfigure(0, weight=1)

# Entry Boxes
log_const_min = ttk.Entry(log_frame)
log_const_min.grid(row=1, padx=5, sticky="ew")
log_const_max = ttk.Entry(log_frame)
log_const_max.grid(row=3, padx=5, sticky="ew")

# Labels
log_const_min_label = ttk.Label(log_frame, text="Min Log Const.").grid(row=0, sticky="ew", padx=5)
log_const_max_label = ttk.Label(log_frame, text="Max Log Const.").grid(row=2, sticky="ew", padx=5)


#Arithmetic Settings
arith_frame = tk.Frame(app, width=200, height=100)
arith_frame.grid_propagate(False)
arith_frame.columnconfigure(0, weight=1)

# Entry Boxes
arith_const_min = ttk.Entry(arith_frame)
arith_const_min.grid(row=1, padx=5, sticky="ew")
arith_const_max = ttk.Entry(arith_frame)
arith_const_max.grid(row=3, padx=5, sticky="ew")

# Labels
arith_const_min_label = ttk.Label(arith_frame, text="Min Arithmetic Const.").grid(row=0, sticky="ew", padx=5)
arith_const_max_label = ttk.Label(arith_frame, text="Max Arithmetic Const.").grid(row=2, sticky="ew", padx=5)

# Geometric Settings
geom_frame = tk.Frame(app, width=200, height=100)
geom_frame.grid_propagate(False)
geom_frame.columnconfigure(0, weight=1)

# Entry Boxes
geom_const_min = ttk.Entry(geom_frame)
geom_const_min.grid(row=1, padx=5, sticky="ew")
geom_const_max = ttk.Entry(geom_frame)
geom_const_max.grid(row=3, padx=5, sticky="ew")

# Labels
geom_const_min_label = ttk.Label(geom_frame, text="Min Geometric Const.").grid(row=0, sticky="ew", padx=5)
geom_const_max_label = ttk.Label(geom_frame, text="Max Geometric Const.").grid(row=2, sticky="ew", padx=5)

# Periodic Settings
period_frame = tk.Frame(app, width=200, height=100)
period_frame.grid_propagate(False)
period_frame.columnconfigure(0, weight=1)

# Entry Boxes
period_const_min = ttk.Entry(period_frame)
period_const_min.grid(row=1, padx=5, sticky="ew")
period_const_max = ttk.Entry(period_frame)
period_const_max.grid(row=3, padx=5, sticky="ew")

# Labels
period_const_min_label = ttk.Label(period_frame, text="Min Periodic Const.").grid(row=0, sticky="ew", padx=5)
period_const_max_label = ttk.Label(period_frame, text="Max Periodic Const.").grid(row=2, sticky="ew", padx=5)

# Sim Settings
frame_sim_sets = tk.Frame(app, width=200, height=200)
frame_sim_sets.grid(column=0, row=1, columnspan=2, sticky="nsew")

frame_sim_sets.columnconfigure(list(range(2)), weight=1, uniform="Silent_Creme")



start_cadence_label = tk.Label(frame_sim_sets, text="Start Cadence (d)").grid(row=0, sticky="ew", padx=5)
start_cadence = ttk.Entry(frame_sim_sets, width=10)
start_cadence.grid(column=0, row=1, sticky="ew", padx=5)

start_offset_label = tk.Label(frame_sim_sets, text="Start Offset (d)").grid(row=2, sticky="ew", padx=5)
start_offset = ttk.Entry(frame_sim_sets, width=10)
start_offset.grid(column=0, row=3, sticky="ew", padx=5)

max_gap_label = tk.Label(frame_sim_sets, text="Max Gap (d)").grid(row=4, sticky="ew", padx=5)
max_gap = ttk.Entry(frame_sim_sets, width=10)
max_gap.grid(column=0, row=5, sticky="ew", padx=5)

num_iterations_label = tk.Label(frame_sim_sets, text="Number of Iterations").grid(row=6, sticky="ew", padx=5)
num_iterations = ttk.Entry(frame_sim_sets, width=10)
num_iterations.grid(column=0, row=7, sticky="ew", padx=5)

#toa file
toa_file_name_label = tk.Label(frame_sim_sets, text="TOA File Name").grid(row=0, column=1, sticky="ew", padx=5)
toa_file_name = ttk.Entry(frame_sim_sets, width=10)
toa_file_name.grid(column=1, row=1, sticky="ew", padx=5)

file_name_label = tk.Label(frame_sim_sets, text="PNG Save File Name").grid(row=2, column=1, sticky="ew", padx=5)
file_name = ttk.Entry(frame_sim_sets, width=10)
file_name.grid(column=1, row=3, sticky="ew", padx=5)


sim_button = ttk.Button(frame_sim_sets, text="Simulate", command=simulate)
sim_button.grid(column=1, row=7, sticky="nsew")

simprogress = ttk.Progressbar(app, orient=tk.HORIZONTAL)
simprogress.grid(column=0, row=3, columnspan=2, sticky="esw")

app.mainloop()
