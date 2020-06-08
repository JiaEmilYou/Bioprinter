#!/usr/bin/python
# -*- coding: utf-8 -*-
# fiting with Python
# if you are using Spyder: press F6 and select  
# "Execute in dedicated Python console" and 
# "Interact with Python console after execution"

'''this program reads data from a text file of pressure and weight versus time and calculates the mean flow rate at each pressure plateau.
Use the output text file of the printing as input
The output includes 3 plots and a text file.
figure1: pressure and flow rate over time
figure2: pressure and mass over time
figure3: mean flow rate at each pressure plateau, if it's a Newtonian fluid, fit a linear relation to it, the slope of the line can be used to calculated the viscosity of the material through a hagen poiseuille equation
text: mean flow rate at each pressure plateau(mean measured pressure)

'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, leastsq
from tkinter import Tk
from tkinter import filedialog
#from tkFileDialog import askopenfilename
import csv
import os


'''this is what you might need to change'''

density=1.0 #density of the bio ink in units of g/mL

#if it is a Newtonian Fluid
IsNewtonianFluid=False #True of False, if True, do a linear fit




'''common problem handling'''

thf=4 #increase this if one pressure plateau gives move than one data point.


'''don't change the following unless you are sure about what you are doing'''



plt.close('all')
# initialization of user interface
root = Tk()
root.withdraw() # we don't want a full GUI, so keep the root window from appearing



# ----------define the fit functions----------
def fitfunc(x, p0):  # for curve_fit
    return p0 * x



#----------read data from a csv (text) file-----------------------------------
path_filename = filedialog.askopenfilename(title="select the data file",filetypes=[("txt file",'*.txt')]) # show an "Open" dialog box and return the path to the selected file
filename = os.path.basename(path_filename)
filename_base, file_extension = os.path.splitext(filename)
path = os.path.dirname(path_filename)
os.chdir(path)

t=[]
p=[]
w=[]

with open(path_filename) as csvfile:
    rowreader = csv.reader(csvfile, delimiter='\t') # please check the delimiter of your data file!!!!!!!!!
    for i, line in enumerate(rowreader):
        print(', '.join(line))
        t.append(np.float(line[0]))
        p.append(np.float(line[1]))
        w.append(np.float(line[2]))
t=np.array(t)
p=np.array(p)
w=np.array(w)


############Calculation, fit and plot########################################
plt.ioff()  #interactive mode "on" means that figures and plots remain open after program has ended
            #interactive mode "off" means that figures and plots need to be closed for the programm to continue

#----------general fonts for plots and figures----------
font = {'family' : 'sans-serif',
        'sans-serif':['Arial'],
        'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=18)
#seaborn.set(font_scale=1.5)





#----------create a new figure window for the flow and pressure versus time ----------
fig1=plt.figure(1, (10, 5))
border_width = 0.2
ax_size = [0+border_width, 0+border_width, 
           1-2*border_width, 1-2*border_width]
ax_a = fig1.add_axes(ax_size)

ax_a.plot(t,p*100, '-', color = '#1f77b4', linewidth=2, label="pressure", zorder=1) #plot the  data
ax_a.set_ylabel('pressure (kPa)', color ='#1f77b4' )
ax_a.tick_params('y', colors='#1f77b4')
plt.ylim(0,np.max(p*100)*1.2)
plt.legend(loc=2, numpoints=1)

ax_b = ax_a.twinx()  
flow=(w[1:]-w[0:-1])/(t[1:]-t[0:-1])     
ax_b.plot(t[0:-1],flow, '.-', color = '#d62728', linewidth=1, markerfacecolor='#d62728', markersize=6.0,markeredgewidth=0, label="flow", zorder=1) #plot the  data
ax_b.set_ylabel('flow (µl/s)',color = '#d62728')
ax_b.tick_params('y', colors='#d62728')    
plt.ylim(0,np.max(flow)*1.15)
plt.legend(loc=1, numpoints=1)

ax_a.set_xlabel('time (s)')
plt.show()
plt.savefig(filename_base + ' flow and pressure versus time.png', dpi=300, format='png')





#----------create a new figure window for the weight and pressure versus time ----------
fig2=plt.figure(2, (10, 6))
border_width = 0.2
ax_size = [0+border_width, 0+border_width, 
           1-2*border_width, 1-2*border_width]
ax_l = fig2.add_axes(ax_size)

ax_l.plot(t,p*100, '-', color = '#1f77b4', linewidth=2, label="pressure", zorder=2) #plot the  data
ax_l.set_ylabel('pressure (kPa)', color ='#1f77b4' )
ax_l.tick_params('y', colors='#1f77b4')
plt.ylim(0,np.max(p*100)*1.2)
plt.legend(loc=2, numpoints=1)

ax_r = ax_l.twinx()
ax_r.plot(t,(w-w[0]), '.-', color = '#d62728', markerfacecolor='#d62728', markersize=5.0,markeredgewidth=0, label="weight", zorder=3) #plot the  data
ax_r.set_ylabel('weight (mg)',color = '#d62728')
ax_r.tick_params('y', colors='#d62728')    
plt.ylim(0,np.max(w-w[0])*1.15)
plt.legend(loc=1, numpoints=1)
             
ax_l.set_xlabel('time (s)')
plt.show()
plt.savefig(filename_base + ' weight and pressure versus time.png', dpi=300, format='png')




#mean flow rate for pressure plateau
dp = p[1:]-p[0:-1]
ix =np.where(p==np.max(p))[0] #find the location of the maximum pressure
th = np.std(dp[0:ix[0]])*thf #estimate the pressure jump between pressure plateaus
ix=np.where(np.abs(dp)>th)[0] #find the time points where the pressure has jumped
i_start = 0
p_mean=[]
flow_mean=[]
for i in ix:
    if i-i_start>2: #at least 3 measurement values at each pressure level are needed for a stable reading
        p_mean.append(100*np.mean(p[i_start:i]))
        flow_mean.append((w[i]-w[i_start])/(t[i]-t[i_start])/density)
        #flow_mean.append((w[i]-w[i_start])/(t[i]-t[i_start])/1.26) #considering the density of glycerol
        plt.figure(1)
        ax_a.fill_between([t[i_start],t[i]], [0,0], [np.max(p*100)*1.2,np.max(p*100)*1.2], facecolor='#d8d8d8', interpolate=True, zorder=0)
    i_start = i+1
p_mean=np.array(p_mean)
flow_mean = np.array(flow_mean)






#----------create a new figure window for the mean flow and pressure  ----------
fig3=plt.figure(3, (4, 4))
border_width = 0.2
ax_size = [0+border_width, 0+border_width,
           1-2*border_width, 1-2*border_width]
ax3 = fig3.add_axes(ax_size)
plt.figure(3)
plt.plot(p_mean,flow_mean,'o', label="data", zorder=1)
plt.xlabel('pressure (kPa)')
plt.ylabel('flow (µl/s)')
plt.show()
#figure 3 continues in the next block if using Newtonian Fluid



#----------Linear fit of Newtonian Fluid and add the fit to figure3----------
if IsNewtonianFluid==True:
    pstart=(10) #initial guess
    pfit, pcov = curve_fit(fitfunc, p_mean, flow_mean, pstart, bounds = ([0],[1000])) #do the fitting
    err = (np.diag(pcov))**0.5 #estimate 1 standard error of the fit parameters
    print("Fit Parameter: p0=%.3f µl/s/kPa+- %.3f " %(pfit[0],err[0]))
    plt.plot(p_mean,fitfunc(p_mean,pfit),'--', color='#d62728', label="lin fit", zorder=0)
    plt.legend(loc=0, numpoints=1)
    plt.text(50,1,'slope = '+"{:.3f}".format(pfit[0])+'µl/s/kPa', fontsize=12)
    plt.show()
plt.savefig(filename_base + ' mean flow versus mean pressure.png', dpi=300, format='png')

# save the mean flow versus mean pressure data to a file
f = open(filename_base + ' mean.txt','w')
for i in range(0,len(p_mean)): 
    f.write(str(p_mean[i])+'\t'+str(flow_mean[i])+'\n')
f.close()
