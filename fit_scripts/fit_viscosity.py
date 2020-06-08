#!/usr/bin/python
# ============================================================================
#            IMPORT GENERAL PACKAGES
# ============================================================================

#check 1. needle parameter  2. pressure used   3. vis_init   4. plot and save


import sys
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, leastsq,least_squares
from tkinter import Tk
from tkinter import filedialog
#from tkFileDialog import askopenfilename
import csv

plt.close('all')
# initialization of user interface
root = Tk()
root.withdraw() # we don't want a full GUI, so keep the root window from appearing



'''this is what you might need to change'''
# ============================================================================
#            IMPORT CUSTOM CLASSES
# ============================================================================
sys.path.insert(0, '/home/emil/Documents/JiaMasterProject/printingscripts/modeltest/classes/') #change the path of your classes folder
from interpolation import Interpolation, Analytical_Viscosity #don't worry if it says 'unresolved reference'
from fluidprofile import Printing_Parameters, Profiles #don't worry if it says 'unresolved reference'




# ============================================================================
#            SET PARAMETERS
# ============================================================================



rchannel = 280/(10**6)       # radius of the cylindrical channel (in m) (inner diameter of the needle)
lchannel = 0.028             #length of the needle (in m)
vis_init=[2.169,57.92,0.729] #vis_0, gamma_c, alpha #initial values for least square iteration








'''don't change the following unless you are sure about what you are doing'''





# --- interpolation parameters ---
gamma0 = 1.0e-6               # gamma0 is the start shear rate for the interpolation
gammaN = 1.0e8                # gammaN is the end shear rate for the interpolation
Ninterpol = 100               # Ninterpol is the number of interpolation intervals [10-1000]
# optional interpolation parameters
gammaStart = gamma0 * 1.0e-2  # gammaStart: Only for plotting - the lower limit of the shear rate
gammaEnd   = gammaN * 1.0e4   # gammaEnd:   Only for plotting - the upper limit of the shear rate
samples    = 10000            # samples:    Only for plotting - the number of calculated points [10000-1000000]



############Plotting and all that stuff########################################
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

#----------read data from a csv (text) file-----------------------------------
filename = []
filename = filedialog.askopenfilename(title="select the data file",filetypes=[("txt file",'*.txt')]) # show an "Open" dialog box and return the path to the selected file

if filename == '':
    print('empty')
    sys.exit()
else:
    pathname=os.path.dirname(filename)
    path_filename = os.path.join(pathname, filename)
    print(path_filename)
    p=[]
    f=[]
    os.chdir(pathname)
#----------read text file------------------------------------    
    with open(filename) as csvfile:
        rowreader = csv.reader(csvfile, delimiter='\t') # please check the delimiter of your data file!!!!!!!!!
        for i, line in enumerate(rowreader):
            print(', '.join(line))
            p.append(np.float(line[0]))
            f.append(np.float(line[1]))
        

p=np.array(p)
f=np.array(f)
f[0]=0
fl_fit=np.zeros(len(p))

def lsqfunc(para,pr,fl): # for leastsq
    # --- analytical model parameters (Carreau-Yasuda model) ---
    eta0   = para[0]          # viscosity in the limit of zero shear rates (in Pa*s)
    etainf = 0.00001          # viscosity in the limit of infinite shear rates (in Pa*s)
    K      = 1/para[1]         # consistency parameter (sometimes "corner shear rate") (in s)
    a1     = para[2]            # first exponent
    a2     = para[2]            # second exponent
    global fl_fit
    fl_fit=np.zeros(len(pr))
    print('eta = %5.2f  gamma_c= %5.2f   a=%5.3f' % (para[0],para[1],para[2]))
    for i,pressure in enumerate(pr):
        pressure = - pressure*1000
        analytical = Analytical_Viscosity(eta0=eta0, etainf=etainf, K=K, a1=a1, a2=a2)
        interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples, gammaStart=gammaStart, gammaEnd=gammaEnd)
        interpol.calculate_interpolation()
        printparams = Printing_Parameters(pressureDifference=pressure, channelLength=lchannel, channelRadius=rchannel)      # with pressure difference and channel length
        fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
        fluidprofiles.calculate_profiles(loud=1)        
        fl_fit[i] = fluidprofiles.flowrate*10**9
        print('pressure %8.1f kPa   flow_meas %6.3f  flow_fit %6.3f' % (-pressure/1000,fl[i],fl_fit[i]))
    print('\n')
    return np.log10(fl+1)-np.log10(fl_fit+1)#minimize the sum of this value at different pressure



vis=least_squares(lsqfunc,vis_init,args=(p,f),bounds = ([0.0001,0.001,0.5],[10000,10**10,1]),ftol=1e-8)
print(vis.x)




#----------create a new figure window----------
fig1=plt.figure(1, (6, 6))
border_width = 0.2
ax_size = [0+border_width+0.05, 0+border_width, 
           1-2*border_width+0.05, 1-2*border_width]
ax1 = fig1.add_axes(ax_size)


#--------do the actual plotting--------------
plt.plot(p,f, 'o', markerfacecolor='#1f77b4', markersize=6.0,markeredgewidth=0, label="data", zorder=3) #plot the  data
#plt.plot(p,fl_fit,'--',color = '#1f77b4',linewidth=2.0,label="fit", zorder=1)
plt.plot(p,fl_fit,'--',color = '#ff7f0e',linewidth=1.5,label="fit", zorder=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('pressure (kPa)')
plt.ylabel('flow (µl/s)')
#plt.xticks([50,100,150,200,250,300,350,400],['50','100','','200','','','','400'])
plt.minorticks_off()
plt.legend(loc=0, numpoints=1)
s = "$\eta_0$ = %0.3f Pa s\n$\gamma_c$ = %0.2f 1/s\n$\\alpha$  = %0.3f" % (vis.x[0], vis.x[1], vis.x[2]) #convert to string
plt.text(80,10,s,fontsize=12)
plt.show()


