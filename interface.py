from tkinter import *
from tkinter import ttk
import threading

import simulating_cadence

app = Tk()
app.title("Pulsar Timer")

def simulate_thread():

    const_args = (start_cadence.get(), start_offset.get(), max_gap.get())
    sim_args = (log_const_min.get(), log_const_max.get(), num_iterations.get())
    print("simulating with following args:")
    print("const_args: ")
    print(const_args)
    print("sim_args: ")
    print(sim_args)

    simulating_cadence.simulate("master_toas.tim", SEQUENCE_TYPE.get(), const_args, sim_args)


def schedule_check(thread):
    if thread.is_alive():
        app.after(100, schedule_check, thread)
    else:
        print("simulation complete")

def forget_all():
    log_frame.grid_forget()
    arith_frame.grid_forget()
    geom_frame.grid_forget()
    period_frame.grid_forget()
    default_frame.grid_forget()

def log_frame_command():
    forget_all()
    log_frame.grid(column=1, row=0)

def arith_frame_command():
    forget_all()
    arith_frame.grid(column=1, row=0)

def geom_frame_command():
    forget_all()
    geom_frame.grid(column=1, row=0)

def period_frame_command():
    forget_all()
    period_frame.grid(column=1, row=0)

def simulate():

    sim_thread = threading.Thread(target=simulate_thread)
    sim_thread.start()
    schedule_check(sim_thread)

frame_radio = Frame(app, width=200, height=200)
frame_radio.grid(column=0, row=0)

default_frame = Frame(app, width=200, height=200)
default_frame.grid(column=1, row=0)



# Create 4 radio buttons in the left frame
SEQUENCE_TYPE = StringVar()
SEQUENCE_TYPE.set('none')
ttk.Radiobutton(frame_radio, text='Logarithmic', variable=SEQUENCE_TYPE, value='logarithmic',command=log_frame_command).grid(column=0, row=0)
ttk.Radiobutton(frame_radio, text='Arithmetic', variable=SEQUENCE_TYPE, value='arithmetic',command=arith_frame_command).grid(column=0, row=1)
ttk.Radiobutton(frame_radio, text='Geometric', variable=SEQUENCE_TYPE, value='geometric',command=geom_frame_command).grid(column=0, row=2)
ttk.Radiobutton(frame_radio, text='Periodic', variable=SEQUENCE_TYPE, value='periodic',command=period_frame_command).grid(column=0, row=3)

# Log Settings
log_frame = Frame(app, width=200, height=200)

# Entry Boxes
log_const_min = ttk.Entry(log_frame, width=10)
log_const_min.grid(column=0, row=1)
log_const_max = ttk.Entry(log_frame, width=10)
log_const_max.grid(column=0, row=3)

# Labels
log_const_min_label = ttk.Label(log_frame, text="Min Log Const.").grid(column=0, row=0)
log_const_max_label = ttk.Label(log_frame, text="Max Log Const.").grid(column=0, row=2)



arith_frame = Frame(app, width=200, height=200)
geom_frame = Frame(app, width=200, height=200)
period_frame = Frame(app, width=200, height=200)



frame_sim_sets = Frame(app, width=200, height=200)
frame_sim_sets.grid(column=0, row=1, columnspan=2)


start_cadence_label = Label(frame_sim_sets, text="Start Cadence (d)").grid(column=0, row=0)
start_cadence = ttk.Entry(frame_sim_sets, width=10)
start_cadence.grid(column=0, row=1)

start_offset_label = Label(frame_sim_sets, text="Start Offset (d)").grid(column=0, row=4)
start_offset = ttk.Entry(frame_sim_sets, width=10)
start_offset.grid(column=0, row=5)

max_gap_label = Label(frame_sim_sets, text="Max Gap (d)").grid(column=0, row=2)
max_gap = ttk.Entry(frame_sim_sets, width=10)
max_gap.grid(column=0, row=3)

num_iterations_label = Label(frame_sim_sets, text="Number of Iterations").grid(column=0, row=6)
num_iterations = ttk.Entry(frame_sim_sets, width=10)
num_iterations.grid(column=0, row=7)

file_name_label = Label(frame_sim_sets, text="File Name").grid(column=1, row=0)
file_name = ttk.Entry(frame_sim_sets, width=10)
file_name.grid(column=1, row=1)


sim_button = ttk.Button(frame_sim_sets, text="Simulate", command=simulate).grid(column=1, row=5)

simprogress = ttk.Progressbar(app, orient=HORIZONTAL).grid(column=0, row=3, columnspan=2)

app.mainloop()