#!/usr/bin/python
# ============================================================================
#            IMPORT GENERAL PACKAGES
# ============================================================================
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
import matplotlib
matplotlib.rc('legend', fontsize=15)

plt.close('all')
# initialization of user interface
root = Tk()
root.withdraw() # we don't want a full GUI, so keep the root window from appearing



'''this is what you might need to change'''
# ============================================================================
#            IMPORT CUSTOM CLASSES
# ============================================================================
sys.path.insert(0, '/home/emil/Documents/JiaMasterProject/printingscripts/modeltest/classes/')
from interpolation import Interpolation, Analytical_Viscosity
from fluidprofile import Printing_Parameters, Profiles

#needle parameters
#parameter of a needle with 230µm diameter
rchannel230 = 115/(10**6)       # radius of the cylindrical channel (in m)
lchannel230 = 0.018 #in m
#parameter of a needle with 560µm diameter
rchannel560 = 280/(10**6)       # radius of the cylindrical channel (in m)
lchannel560 = 0.028 #in m

vis_init=[2.169,57.92,0.729] #vis_0, gamma_c, alpha  #initial values for least square iteration








'''don't change the following unless you are sure about what you are doing'''







############Plotting and all that stuff########################################
plt.ioff()  #interactive mode "on" means that figures and plots remain open after program has ended
            #interactive mode "off" means that figures and plots need to be closed for the programm to continue


#----------read data from a csv (text) file-----------------------------------
def readflowpressure(title):
    filename = []
    filename = filedialog.askopenfilename(title="select the data file"+title,filetypes=[("txt file",'*.txt')]) # show an "Open" dialog box and return the path to the selected file
    
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
            
    #            if np.float(line[1])>0.5: #only accept data with a flow rate of larger than 0.5 ul/s
    #                p.append(np.float(line[0]))
    #                f.append(np.float(line[1]))
    p=np.array(p)
    f=np.array(f)
    f[0]=0
    fl_fit=np.zeros(len(p))
    return(np.vstack((p,f,fl_fit)))

data560=readflowpressure(" 560 needle")
data230=readflowpressure(" 230 needle")

# ============================================================================
#            SET PARAMETERS
# ============================================================================
# --- interpolation parameters ---
gamma0 = 1.0e-6               # gamma0 is the start shear rate for the interpolation
gammaN = 1.0e8                # gammaN is the end shear rate for the interpolation
Ninterpol = 150               # Ninterpol is the number of interpolation intervals [10-1000]
# optional interpolation parameters
gammaStart = gamma0 * 1.0e-2  # gammaStart: Only for plotting - the lower limit of the shear rate
gammaEnd   = gammaN * 1.0e4   # gammaEnd:   Only for plotting - the upper limit of the shear rate
samples    = 10000            # samples:    Only for plotting - the number of calculated points [10000-1000000]


def lsqfunc(para,prg,flg,prs,fls): # for leastsq
    # --- analytical model parameters (Carreau-Yasuda model) ---
    eta0   = para[0]          # viscosity in the limit of zero shear rates (in Pa*s)
    #eta0   = 8          # viscosity in the limit of zero shear rates (in Pa*s)
    etainf = 0.00001          # viscosity in the limit of infinite shear rates (in Pa*s)
    K      = 1/para[1]         # consistency parameter (sometimes "corner shear rate") (in s)
    a1     = para[2]            # first exponent
    a2     = para[2]            # second exponent

    global fl_fit
    fl=np.concatenate((data560[1,],data230[1,]))
    fl_fit=np.zeros(len(prg)+len(prs))
    print('eta=%5.2f  gamma_c= %5.2f   a=%5.3f' % (para[0],para[1],para[2]))
    #print('gamma_c= %5.2f   a=%5.3f' % (para[0],para[1]))
    for i,pressure in enumerate(prg):
        pressure = - pressure*1000
        analytical = Analytical_Viscosity(eta0=eta0, etainf=etainf, K=K, a1=a1, a2=a2)
        interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples, gammaStart=gammaStart, gammaEnd=gammaEnd)
        interpol.calculate_interpolation()
        printparams = Printing_Parameters(pressureDifference=pressure, channelLength=lchannel560, channelRadius=rchannel560)      # with pressure difference and channel length
        fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
        fluidprofiles.calculate_profiles(loud=1)        
        fl_fit[i] = fluidprofiles.flowrate*10**9
        print('pressure %8.1f kPa   flow_meas %6.1f  flow_fit %6.1f' % (-pressure/1000,flg[i],fl_fit[i]))
    for i,pressure in enumerate(prs):
        pressure = - pressure*1000
        analytical = Analytical_Viscosity(eta0=eta0, etainf=etainf, K=K, a1=a1, a2=a2)
        interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples, gammaStart=gammaStart, gammaEnd=gammaEnd)
        interpol.calculate_interpolation()
        printparams = Printing_Parameters(pressureDifference=pressure, channelLength=lchannel230, channelRadius=rchannel230)      # with pressure difference and channel length
        fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
        fluidprofiles.calculate_profiles(loud=1)        
        fl_fit[len(prg)+i] = fluidprofiles.flowrate*10**9
        print('pressure %8.1f kPa   flow_meas %6.1f  flow_fit %6.1f' % (-pressure/1000,fls[i],fl_fit[len(prg)+i]))
    print('\n')
    return np.log10(fl+0.001)-np.log10(fl_fit+0.001)#minimize the sum of this value at different pressure





vis=least_squares(lsqfunc,vis_init,args=(data560[0,],data560[1,],data230[0],data230[1,]),bounds = ([0.1,0.001,0.5],[100,10**10,1]),ftol=1e-8)
#vis=least_squares(lsqfunc,vis_init,args=(data560[0,],data560[1,],data230[0],data230[1,]),bounds = ([0.001,0.5],[10**10,1]),ftol=0.000001)
print(vis.x)


#----------create a new figure window----------



pres_fit_560 = [data560[0,1]*0.8]
while pres_fit_560[-1]<data560[0,-1]*1.1:
    pres_fit_560.append(pres_fit_560[-1]*1.1)
pres_fit_560 = np.asarray(pres_fit_560)
flow_fit_560 =  0*pres_fit_560

pres_fit_230 = [data230[0,1]*0.8]
while pres_fit_230[-1]<data230[0,-1]*1.1:
    pres_fit_230.append(pres_fit_230[-1]*1.1)
pres_fit_230 = np.asarray(pres_fit_230)
flow_fit_230 =  0*pres_fit_230


for i,p in enumerate(pres_fit_560):
        p = - p*1000
        analytical = Analytical_Viscosity(eta0=vis.x[0], etainf=0.00001, K=1/vis.x[1], a1=vis.x[2], a2=vis.x[2])
        interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples, gammaStart=gammaStart, gammaEnd=gammaEnd)
        interpol.calculate_interpolation()
        printparams = Printing_Parameters(pressureDifference=p, channelLength=lchannel560, channelRadius=rchannel560)      # with pressure difference and channel length
        fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
        fluidprofiles.calculate_profiles(loud=1)        
        flow_fit_560[i] = fluidprofiles.flowrate*10**9
        
for i,p in enumerate(pres_fit_230):
        p = - p*1000
        analytical = Analytical_Viscosity(eta0=vis.x[0], etainf=0.00001, K=1/vis.x[1], a1=vis.x[2], a2=vis.x[2])
        interpol = Interpolation(gamma0=gamma0, gammaN=gammaN, Ninterpol=Ninterpol, analytical=analytical, samples=samples, gammaStart=gammaStart, gammaEnd=gammaEnd)
        interpol.calculate_interpolation()
        printparams = Printing_Parameters(pressureDifference=p, channelLength=lchannel230, channelRadius=rchannel230)      # with pressure difference and channel length
        fluidprofiles = Profiles(interpolation=interpol, printingParameters=printparams, samples=samples)
        fluidprofiles.calculate_profiles(loud=1)        
        flow_fit_230[i] = fluidprofiles.flowrate*10**9
        
fig1=plt.figure(1, (7, 7))
border_width = 0.2
ax_size = [0+border_width+0.05, 0+border_width, 
           1-2*border_width+0.05, 1-2*border_width]
ax1 = fig1.add_axes(ax_size)
#--------do the actual plotting--------------
plt.plot(data560[0,],data560[1,], 'o', markerfacecolor='#2ca02c', markersize=7.0,markeredgewidth=0, label="560µm data", zorder=3) #plot the  data
#plt.plot(p,fl_fit,'--',color = '#1f77b4',linewidth=2.0,label="fit", zorder=1)
plt.plot(pres_fit_560,flow_fit_560,'--',color = '#ff7f0e',linewidth=1.5, zorder=1)
         
plt.plot(data230[0,],data230[1,], 'o', markerfacecolor='#1f77b4', markersize=7.0,markeredgewidth=0, label="230µm data", zorder=3) #plot the  data
#plt.plot(p,fl_fit,'--',color = '#1f77b4',linewidth=2.0,label="fit", zorder=1)
plt.plot(pres_fit_230,flow_fit_230,'--',color = '#ff7f0e',linewidth=1.5,label="fit", zorder=1)
         

s = "$\eta_0$ = %0.3f Pa s\n$\gamma_c$ = %0.2f 1/s\n$\\alpha$  = %0.3f" % (vis.x[0], vis.x[1], vis.x[2]) #convert to string
#plt.text(100,0.2,s,fontsize=18)
plt.text(12,0.02,s,fontsize=18) #2.5
#plt.savefig(r'flow_vs_pressure_global_fit.png', dpi = 300, format='png')
#plt.show()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('pressure (kPa)',size=24)
plt.ylabel('flow (µl/s)',size=24)
plt.tick_params(labelsize=20)
#plt.xticks([10,15,20,30,40,60,80,120,160,240,320],['10','','20','','40','','80','','160','','320'])
#plt.ylim([0.05,250])
plt.ylim([0.01,150])
plt.minorticks_off()
plt.legend(loc=2, numpoints=1)
plt.show()

    
    
