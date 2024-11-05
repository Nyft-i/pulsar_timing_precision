# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:17:16 2024

@author: laure
"""

import numpy as np 
import matplotlib.pyplot as plt

#results for pulsars 1 to 4 from manual tempo2 fitting
#order goes log > geo > periodic
pulsar_one_off = np.array([[4.70798e-7,0.00000000000138144530,-2.099997e-15, 2.0778145875309407317e-20],
                           [4.7137972923489013029e-07,0.00000000000140384872,-2.1000108110149476731e-15,2.1268451147202375242e-20],
                           [4.7000861410119862005e-07,0.00000000000139039028,-2.0999882135220762017e-15,2.1118070273141016292e-20]])
                  
pulsar_one = [4.7e-07,-2.1e-15]

pulsar_two_off = np.array([[8.9782936805054591118e-08,1.6018365554312269025e-13,-3.5000014157064762785e-15, 2.0529318766903285723e-21],
                           [8.8361679472061699872e-08,1.5665785593660629774e-13,-3.5000019831927997162e-15,2.0253192268613766227e-21],
                           [8.8028808169033183034e-08,1.5578043783886959696e-13,-3.4999991848534033674e-15,2.0124280403121040857e-21]])
                  
pulsar_two = [8.88e-08,-3.5e-15]

pulsar_three_off = np.array([[5.6547269534285927437e-07,0.00000000000315605664,-3.1399664266069330584e-15,3.4378211330812529452e-20],
                             [5.6353505421336287701e-07,0.00000000000292105549,-3.1399931669298158158e-15,3.2086408932266411019e-20],
                             [5.653986464061750689e-07,0.00000000000296964249,-3.1399048331740586785e-15,3.2464137359067716059e-20]])
                  
pulsar_three = [5.65e-7,-3.14e-15]

pulsar_four_off = np.array([[2.3648569933575638997e-07,3.1268763052173749551e-13,-1.1899982723159473428e-15,5.5309339391854150342e-21],
                            [2.3682041389808351409e-07,3.2244024112191226885e-13,-1.1900078623114164471e-15,5.7065638810490373625e-21],
                            [2.3712337684630463296e-07,3.2159399252618939285e-13,-1.1899994549210307983e-15,5.6977393136013276436e-21]])
                  
pulsar_four = [2.37e-7,-1.19e-15]

fig = plt.figure(figsize=(16, 4))
gs = fig.add_gridspec(2, 2, wspace = 0, hspace = 0)
axs = gs.subplots(sharey=True)

fig.suptitle("")


axs[0].errorbar(pulsar_one_off[0,0]-pulsar_one[0], pulsar_one_off[0,2]-pulsar_one[1], xerr=pulsar_one_off[0,1], yerr=pulsar_one_off[0,3], fmt='x', label="logarithmic")
axs[0].errorbar(pulsar_one_off[1,0]-pulsar_one[0], pulsar_one_off[1,2]-pulsar_one[1], xerr=pulsar_one_off[1,1], yerr=pulsar_one_off[1,3], fmt='x', label="geometric")
axs[0].errorbar(pulsar_one_off[2,0]-pulsar_one[0], pulsar_one_off[2,2]-pulsar_one[1], xerr=pulsar_one_off[2,1], yerr=pulsar_one_off[2,3], fmt='x', label="periodic")

axs[1].errorbar(pulsar_two_off[0,0]-pulsar_two[0], pulsar_two_off[0,2]-pulsar_two[1], xerr=pulsar_two_off[0,1], yerr=pulsar_two_off[0,3], fmt='x', label="logarithmic")
axs[1].errorbar(pulsar_two_off[1,0]-pulsar_two[0], pulsar_two_off[1,2]-pulsar_two[1], xerr=pulsar_two_off[1,1], yerr=pulsar_two_off[1,3], fmt='x', label="geometric")
axs[1].errorbar(pulsar_two_off[2,0]-pulsar_two[0], pulsar_two_off[2,2]-pulsar_two[1], xerr=pulsar_two_off[2,1], yerr=pulsar_two_off[2,3], fmt='x', label="periodic")

axs[2].errorbar(pulsar_three_off[0,0]-pulsar_three[0], pulsar_three_off[0,2]-pulsar_three[1], xerr=pulsar_three_off[0,1], yerr=pulsar_three_off[0,3], fmt='x', label="logarithmic")
axs[2].errorbar(pulsar_three_off[1,0]-pulsar_three[0], pulsar_three_off[1,2]-pulsar_three[1], xerr=pulsar_three_off[1,1], yerr=pulsar_three_off[1,3], fmt='x', label="geometric")
axs[2].errorbar(pulsar_three_off[2,0]-pulsar_three[0], pulsar_three_off[2,2]-pulsar_three[1], xerr=pulsar_three_off[2,1], yerr=pulsar_three_off[2,3], fmt='x', label="periodic")

axs[3].errorbar(pulsar_four_off[0,0]-pulsar_four[0], pulsar_four_off[0,2]-pulsar_four[1], xerr=pulsar_four_off[0,1], yerr=pulsar_four_off[0,3], fmt='x', label="logarithmic")
axs[3].errorbar(pulsar_four_off[1,0]-pulsar_four[0], pulsar_four_off[1,2]-pulsar_four[1], xerr=pulsar_four_off[1,1], yerr=pulsar_four_off[1,3], fmt='x', label="geometric")
axs[3].errorbar(pulsar_four_off[2,0]-pulsar_four[0], pulsar_four_off[2,2]-pulsar_four[1], xerr=pulsar_four_off[2,1], yerr=pulsar_four_off[2,3], fmt='x', label="periodic")


