#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 14:24:43 2019

@author: emil
"""
import sys

sys.path.insert(0, '/home/emil/Documents/JiaMasterProject/modelingscripts/classes/')
from interpolation import Interpolation, Analytical_Viscosity
from fluidprofile import Printing_Parameters, Profiles
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

import pylustrator #pylustrator helps you adjust your image easily, ask richard about pylustrator
pylustrator.start()

font = {'family': 'sans-serif',
        'sans-serif': ['Arial'],
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
matplotlib.rc('legend', fontsize=12)
matplotlib.rc('axes', titlesize=18)


'''this is what you might need to change'''
vis_para=[2.424,42.19,0.677]  # your viscosity parameters from fitting

needle=[]
pressure=[]

#add information for a needle
needle.append([2.8e-2, 2.8e-4])# needle length, needle radius in meter
pressure.append([0.15,0.3,0.45,0.6,0.75,0.9])#printing pressure for this needle

#add information for a needle
needle.append([1.8e-2, 1.15e-4])# needle length, needle radius in meter
pressure.append([1,1.5,2,2.5,3,3.5])#printing pressure for this needle



'''don't change the following unless you are sure about what you are doing'''




# Profiles.print_usage()
# ============================================================================
#            SET PARAMETERS
# ============================================================================
# --- interpolation parameters ---
gamma0 = 1.0e-6  # gamma0 is the start shear rate for the interpolation
gammaN = 1.0e8  # gammaN is the end shear rate for the interpolation
Ninterpol = 150  # Ninterpol is the number of interpolation intervals [10-1000]
# optional interpolation parameters
gammaStart = gamma0 * 1.0e-2  # gammaStart: Only for plotting - the lower limit of the shear rate
gammaEnd = gammaN * 1.0e4  # gammaEnd:   Only for plotting - the upper limit of the shear rate
samples = 10000  # samples:    Only for plotting - the number of calculated points [10000-1000000]

# --- analytical model parameters (Carreau-Yasuda model) ---
# visc = etainf + (eta0 - etainf) / ( (1.0 + (K * gamma)**(a1) )**(a2*1.0/a1) )
eta0 = vis_para[0]  # viscosity in the limit of zero shear rates (in Pa*s)
etainf = 1.0e-7  # viscosity in the limit of infinite shear rates (in Pa*s)
K = 1 / vis_para[1]  # consistency parameter (sometimes "corner shear rate") (in s)
a1 = vis_para[2]  # first exponent
a2 = vis_para[2]  # second exponent

# ============================================================================
#            PERFORM INTERPOLATION AND PROFILE CALCULATION
# ============================================================================ 
# --- Define analytical viscosity model ---
analytical = Analytical_Viscosity(eta0=eta0, etainf=etainf, K=K, a1=a1, a2=a2)

# --- Perform interpolation ---
interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples,
                         gammaStart=gammaStart, gammaEnd=gammaEnd)
interpol.calculate_interpolation()


def calculateprofile_pgrad(p, needle):  # pressure in bar
    pgrad = -p * 1e5 / needle[0]  # pressure 1.8e-2
    printparams = Printing_Parameters(pressureGradient=pgrad, channelRadius=needle[1])
    fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
    fluidprofiles.calculate_profiles(loud=1)
    disc_flowrate = 0
    shear_rate_list = []
    flowrate_i_list = []
    for i in range(1, fluidprofiles.samples + 1, 1):
        ringarea = math.pi * (fluidprofiles.data[0][i] * fluidprofiles.data[0][i] - fluidprofiles.data[0][i - 1] *
                              fluidprofiles.data[0][i - 1])  # area of every ring
        flowrate_i = ringarea * fluidprofiles.data[1][i]  # ringarea*velocity=flow in the ringarea
        disc_flowrate += flowrate_i  # sum of the flow
        shearrate_i = fluidprofiles.data[2][i]
        shear_rate_list.append(shearrate_i)  # list of shear stress in every ring
        flowrate_i_list.append(flowrate_i)  # list of flow rate
    shearratedist = np.ones([2, len(shear_rate_list)])
    shearratedist[0,] = shear_rate_list
    shearratedist[1,] = flowrate_i_list
    dlogshearrate = [1]  # standardize under a log scale
    for i in range(1, len(shearratedist[0,]) - 1):
        dlogshearrate.append(math.log(shearratedist[0, i + 1]) - math.log(shearratedist[0, i]))

    for i in range(1, len(shearratedist[0,]) - 1):
        shearratedist[1, i] /= dlogshearrate[i]
    normsum = sum(shearratedist[1, :])
    # print(normsum)
    for i in range(len(shearratedist[1,])):
        shearratedist[1, i] /= normsum  # the sum of data points is 1 here ######dots sum up to 1
    return (shearratedist)

plotdata=[]
for i in range(len(needle)):
    plotdata.append([])
    for p in pressure[i]:
        plotdata[-1].append(calculateprofile_pgrad(p,needle[i]))


color=['r','g','b','k']

fig = plt.figure(1, (6, 6))
axe = fig.add_axes()
plt.xscale('log')

linestyle = '-'


for i in range(len(needle)):
    for p in range(len(plotdata[i])):
        plt.plot(plotdata[i][p][0], np.cumsum(plotdata[i][p][1]), linestyle, color=color[i], linewidth=2) #'Hello, {}'.format(name)
        plt.text(np.percentile(plotdata[i][p][0],100/(len(needle)+1)*(i+1))+20, np.percentile(np.cumsum(plotdata[i][p][1]),100/(len(needle)+1)*(i+1)+20),'r={}µm p={}bar'.format(round(needle[i][1]*1000000,3),pressure[i][p]), color=color[i],fontsize=10).set_rotation(80)


plt.xlim(1, 100000)
plt.xlabel('Shear rate (1/s)')
plt.ylabel('Cumulative probability')


plt.legend(loc=2, numpoints=1)

plt.show()





