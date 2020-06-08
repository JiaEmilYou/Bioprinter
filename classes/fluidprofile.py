#!/usr/bin/python

# ============================================================================
#            IMPORT PACKAGES
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
import math

from interpolation import Interpolation


# ============================================================================
#            TODO-LIST
# ============================================================================
# - include iterative calculation of pGrad from given flow rate
# - include calculation of flow rate from profile
# - plot from -R:R ? rather mirror plot from R:0:R


# ============================================================================
#            CLASS DEFINITIONS
# ============================================================================
class Printing_Parameters : 
    """PrintingParameters class contains the printing parameters G (pressure gradient, negative) and R (channel radius)"""

    def __init__(self, pressureGradient=None, channelRadius=None, pressureDifference=None, channelLength=None, flowrate=None) :
        """Stores the given parameters. Not all need to be provided"""
        # Check for correct number of arguments
        if channelRadius == None or (pressureGradient == None and (pressureDifference == None or channelLength == None) and flowrate == None) :
            self.print_usage()
        else :
            self.R = channelRadius
            # Either pressure gradient or pressure difference and channellength or flowrate are needed
            if pressureGradient != None :
                self.G        = pressureGradient
                self.flowrate = None
            elif ( pressureDifference != None and channelLength != None) :
                if pressureGradient != None :
                    self.print_usage()
                else :
                    self.G        = pressureDifference * (1.0 / channelLength)
                    self.flowrate = None
            elif flowrate != None :
                if pressureGradient != None :
                    self.print_usage()
                else :                
                    self.flowrate = flowrate
                    self.G        = None
    
    
    @staticmethod
    def print_usage() :
        """Method used to print input instructions."""
        print(" --- Printing_Parameters class usage ---")
        print("\nPrinting_Parameters() class requires at least two parameters:")
        print(" - channelRadius :\t radius of the cylindrical channel (in m)")
        print(" - pressureGradient : \t pressure gradient (negative) along the channel (in Pa/m)")
        print(" - Alternatively to pressureGradient, one of the following parameters must provided:")
        print("    - flowrate :\t\t flow rate in the channel (in m^3/s) (Note: the flow rate is iteratively converted into a pressure gradient, which takes some time (use option 'loud=1' for output))")
        print("    - pressureDifference :\t pressure drop along the channel of length (negative)")
        print("      and channelLength : \t the length of the channel\n\n")

      
class Profiles : 
    """Profiles class contains functions and parameters that define the solution for the shear rate and velocity profile in a channel"""
    
    def __init__(self, interpolation=None, printingParameters=None, samples=10000) :
        """"Initialize the parameters necessary for calculating the semi-analytical solution"""
        if interpolation == None or printingParameters == None :
            self.print_usage()
        else :
            # Interpolated viscosity model
            self.interpol = interpolation
            # Printing parameters    
            self.printparams = printingParameters
            # Save the number of interpolation intervals separately    
            self.Ninterpol = self.interpol.N
            # Initialize the number of relevant interpolation intervals
            self.Nchannel = self.Ninterpol
            # Initialize lists that store the parameters of the analytical solution
            self.r = np.zeros(self.Ninterpol+3)
            self.c = np.zeros(self.Ninterpol+2)
            # Number of intervals necessary to calculate velocity profile
            self.k = 0
            # The plotrange in terms of intervals: 0, ..., k
            self.krange = range(0,self.k+1,1)
            # Matrix to store the calculated values (for quicker printing and plotting)        
            self.samples     = samples
            self.samplerange = range(0,self.samples+1,1)
            self.data        = np.zeros((5, self.samples+1)) # radial position, velocity, shear rate, viscosity, interval (of interpolation)
            # Average values
            self.flowrate      = None
            self.avg_velocity  = None
            self.avg_shearrate = None
            self.avg_viscosity = None
        
    
    @staticmethod
    def print_usage() :
        """Method used to print input instructions."""
        print(" --- Profiles class usage ---")
        print("\nProfiles() class requires at least two parameters:")
        print(" - interpolation :\t\t an instance of the Interpolation() class")
        print(" - printingParameters : \t an instance of the Printing_Parameters() class")
        print(" - (optional) samples :\t\t number of points calculated for plotting and data saving\n")
        print("The following methods must be called:")
        print(" - Profiles() :\t\t\t with all necessary parameters")
        print(" - calculate_profiles() :\t performs the interpolation \n")
        print("The following methods/variables are helpful:")
        print(" - save_profiles(filename) : \t save the profile data to a file (radial position, velocity, shear rate, viscosity)")
        print(" - plot_velocity(options) :\t plot the velocity profile")
        print(" - plot_shearrate(options) :\t plot the shear rate profile")
        print(" - plot_viscosity(options) :\t plot the viscosity profile")
        print("   Options are:")
        print("   - draw_intervals=<int> : \t indicate the interval boundaries with lines")
        print("   - save_figure=<filename> : \t save the plot to a file\n\n")        
 
     
    def find_k(self) :
        """Determine the interval that includes the physical boundary"""
        k = 0
        while ( k <= self.Ninterpol and self.r[k] < self.printparams.R) :
            k += 1
        self.k = k
        self.krange = range(0,self.k+1,1) # 0, ..., k
 
       
    def get_last_interval(self) :
        """Return the last interval that includes the physical boundary (index k)"""
        self.find_k()
        return self.k
        
    def get_max_velocity(self) :
        """Return the maximum velocity"""
        return self.c[0]
        
    def get_max_shearrate(self) :
        """Return the maximum shear rate"""
        return self.calc_shearrate(self.printparams.R, self.k)


    def find_interval(self,r) :
        """Determine the index of the interval that includes the radial position r"""
        i = 0
        while( i < self.Ninterpol and self.r[i] < r) :
            i += 1
        return i
   
     
    def calc_ri(self) :
        """Calculate the interval positions"""
        for i in self.interpol.fullrange :
            self.r[i] = - 2.0 * 1.0/self.printparams.G * self.interpol.K[i] * (self.interpol.gammai[i]**(self.interpol.n[i]))        
    
        
    def calc_ci(self) :
        """Calculate the integrations constants of the velocity field"""
        self.find_k()
        self.c[self.k] = (- 0.5 * self.printparams.G * 1.0/self.interpol.K[self.k])**(1.0/self.interpol.n[self.k]) * self.interpol.n[self.k] * 1.0/(self.interpol.n[self.k] + 1.0) * self.printparams.R**(1.0 + 1.0/self.interpol.n[self.k])
        # The other integration constants are calculated using the already calculated oines
        #for i in range(self.k,-1,-1) : # count backwards from k-1
        for i in range(0,self.k,1) :
            #self.c[i] = self.c[i+1] - self.r[i] * self.interpol.gammai[i] * ( self.interpol.n[i+1]*1.0/(self.interpol.n[i+1]+1.0) - self.interpol.n[i]*1.0/(self.interpol.n[i]+1.0))
            self.c[i] = self.c[self.k]
            for j in range(i,self.k,1) :
                self.c[i] -= self.r[j] * self.interpol.gammai[j] * (self.interpol.n[j+1]/(1.0+self.interpol.n[j+1]) - self.interpol.n[j]/(1.0+self.interpol.n[j])) 
   
                 
    def calculate_profiles(self, loud=None, epsilon=1.0e-10) : 
        """Calculate all coefficients necessary for plotting"""
        if self.printparams.G == None and self.printparams.flowrate != None :
            self.find_pressureGradient(self.printparams.flowrate, loud=loud, epsilon=epsilon)
        self.calc_ri()
        self.calc_ci()
        self.calc_profile_data()
 #       self.flowrate = self.calc_flowrate()
        self.calc_averages()
        
     
     
    def calc_velocity(self,r,i) :
        """Calculate the velocity in interval i at position r"""
        return ( - (- 0.5 * self.printparams.G * 1.0/self.interpol.K[i])**(1.0/self.interpol.n[i]) * self.interpol.n[i] * 1.0/(self.interpol.n[i] + 1.0) * r**(1.0 + 1.0/self.interpol.n[i]) ) + self.c[i]
 
 
    def calc_shearrate(self,r,i) : 
        """Calculate the shear rate in interval i at position r"""
        return (- 0.5 * r * self.printparams.G * 1.0/self.interpol.K[i])**(1.0/self.interpol.n[i])
     
     
    def calc_viscosity(self,r,i) :
        """Calculate the viscosity in interval i at position r"""
        return self.interpol.K[i] * (self.calc_shearrate(r,i))**(self.interpol.n[i]-1.0)
    
    
    def calc_flowrate(self) :
        """Calculates the flow rate for the approximated flow profile"""
        flowrate = 0
        # Add zeroth interval
        if self.k == 0 :
            flowrate += 0.5*self.c[0]*self.printparams.R*self.printparams.R - (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*self.interpol.n[0]*1.0/(1.0+self.interpol.n[0])*1.0/(3.0+1.0/self.interpol.n[0])*(self.printparams.R**(3.0+1.0/self.interpol.n[0]))
        else :
            flowrate += 0.5*self.c[0]*self.r[0]*self.r[0] - (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*self.interpol.n[0]*1.0/(1.0+self.interpol.n[0])*1.0/(3.0+1.0/self.interpol.n[0])*(self.r[0]**(3.0+1.0/self.interpol.n[0]))
        
        # Add k-th interval
        if self.k != 0 :
            flowrate += 0.5*self.c[self.k]*(self.printparams.R**2 - self.r[self.k-1]**2 ) - (-0.5*self.printparams.G*1.0/self.interpol.K[self.k])**(1.0/self.interpol.n[self.k])*self.interpol.n[self.k]*1.0/(1.0+self.interpol.n[self.k])*1.0/(3.0+1.0/self.interpol.n[self.k])*(self.printparams.R**(3.0+1.0/self.interpol.n[self.k]) - self.r[self.k-1]**(3.0+1.0/self.interpol.n[self.k]))
            # Add the intermediate intervals
            for j in range(0,self.k,1) :
                flowrate += 0.5*self.c[j]*(self.r[j]**2 - self.r[j-1]**2 ) - (-0.5*self.printparams.G*1.0/self.interpol.K[j])**(1.0/self.interpol.n[j])*self.interpol.n[j]*1.0/(1.0+self.interpol.n[j])*1.0/(3.0+1.0/self.interpol.n[j])*(self.r[j]**(3.0+1.0/self.interpol.n[j]) - self.r[j-1]**(3.0+1.0/self.interpol.n[j]))
        flowrate *= 2.0 * math.pi
        self.flowrate = flowrate
        return flowrate
        
        
    def find_pressureGradient(self, flowrate, epsilon=None, maxiter=10000, loud=None) :
        """Calculate the pressure gradient from the flow rate iteratively until the relative error is less than epsilon"""
        # Inital values for pressure gradient intervals
        Glow = -1.0e0
        Ghigh = -1.0e1
        # Calculate initial flow rates
        self.printparams.G = Glow
        self.calculate_profiles()
        Qlow = self.calc_flowrate()       
        self.printparams.G = Ghigh
        self.calculate_profiles()
        Qhigh = self.calc_flowrate()
        # Find interval, such that Qlow < flowrate < Qhigh
        while Qlow < flowrate and Qhigh < flowrate :
            if loud != None :
                print("Intervals too low", Glow, Ghigh)
            Glow  = Ghigh
            self.printparams.G = Glow
            self.calculate_profiles()
            Qlow = self.calc_flowrate()           
            Ghigh = 10*Ghigh
            self.printparams.G = Ghigh
            self.calculate_profiles()
            Qhigh = self.calc_flowrate()
        while Qlow > flowrate and Qhigh > flowrate :
            if loud != None :
                print("Intervals too high", Glow, Ghigh)
            Ghigh = Glow
            self.printparams.G = Ghigh
            self.calculate_profiles()
            Qhigh = self.calc_flowrate()        
            Glow  = Glow/10.0
            self.printparams.G = Glow
            self.calculate_profiles()
            Qlow = self.calc_flowrate()           
        # Do bisection until desired accuracy is reached
        Gmid = 0.5 * (Ghigh + Glow)
        self.printparams.G = Gmid
        self.calculate_profiles()
        Qmid = self.calc_flowrate()        
        eps = (Qmid - flowrate)/flowrate
        index = 1
        if loud != None :
            print("iter\tGlow\t\tGhigh\t\tQlow\t\tQmid\t\tQhigh\t\tepsilon")
        while math.fabs(eps) > epsilon :
            Gmid = 0.5 * (Ghigh + Glow)
            self.printparams.G = Gmid
            self.calculate_profiles()
            Qmid = self.calc_flowrate()
            if Qmid > flowrate :
                Ghigh = Gmid
            elif Qmid <= flowrate :
                Glow = Gmid
            self.printparams.G = Glow
            self.calculate_profiles()
            Qlow = self.calc_flowrate()
            self.printparams.G = Ghigh
            self.calculate_profiles()
            Qhigh = self.calc_flowrate()
            eps = (Qmid - flowrate)/flowrate
            if loud != None :
                print("%s\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E"%(index,Glow,Ghigh,Qlow,Qmid,Qhigh,eps))
            index += 1
            if index >= maxiter :
                print("Error: Conversion of flowrate to pressure gradient didnt converge after %s iterations."%(maxiter))
                break
        # Set final pressure gradient between Glow and Ghigh
        Gmid = 0.5 * (Ghigh + Glow)
        self.printparams.G = Gmid
        self.calculate_profiles()
        Qmid = self.calc_flowrate()
        if loud != None :
            print("Conversion of flowrate to pressure gradient converged after %s iterations."%(index))
            print("flowrate (set): %.6E"%(flowrate))
            print("flowrate (get): %.6E"%(Qmid))
        
                
    def calc_averages(self) :
        """Calculate the average velocity, shear rate and viscosity in the channel"""
        area      = math.pi * self.printparams.R * self.printparams.R
        flowrate  = 0
        avg_shear = 0
        avg_visc  = 0
        # Add zeroth interval
        if self.k == 0 :
            flowrate  += 0.5*self.c[0]*self.printparams.R*self.printparams.R - (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*self.interpol.n[0]*1.0/(1.0+self.interpol.n[0])*1.0/(3.0+1.0/self.interpol.n[0])*(self.printparams.R**(3.0+1.0/self.interpol.n[0]))
            avg_shear += (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*1.0/(2.0+1.0/self.interpol.n[0])*(self.printparams.R**(2.0+1.0/self.interpol.n[0]) - self.r[0]**(2.0+1.0/self.interpol.n[0]))
            avg_visc  += self.interpol.K[0] * (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0-1.0/self.interpol.n[0])*1.0/(3.0-1.0/self.interpol.n[0])*(self.printparams.R**(3.0-1.0/self.interpol.n[0]))
        else :
            flowrate += 0.5*self.c[0]*self.r[0]*self.r[0] - (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*self.interpol.n[0]*1.0/(1.0+self.interpol.n[0])*1.0/(3.0+1.0/self.interpol.n[0])*(self.r[0]**(3.0+1.0/self.interpol.n[0]))
            avg_shear += (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0/self.interpol.n[0])*1.0/(2.0+1.0/self.interpol.n[0])*self.r[0]**(2.0+1.0/self.interpol.n[0])
            avg_visc  += self.interpol.K[0] * (-0.5*self.printparams.G*1.0/self.interpol.K[0])**(1.0-1.0/self.interpol.n[0])*1.0/(3.0-1.0/self.interpol.n[0])*(self.r[0]**(3.0-1.0/self.interpol.n[0]))
        # Add k-th interval
        if self.k != 0 :
            flowrate += 0.5*self.c[self.k]*(self.printparams.R**2 - self.r[self.k-1]**2 ) - (-0.5*self.printparams.G*1.0/self.interpol.K[self.k])**(1.0/self.interpol.n[self.k])*self.interpol.n[self.k]*1.0/(1.0+self.interpol.n[self.k])*1.0/(3.0+1.0/self.interpol.n[self.k])*(self.printparams.R**(3.0+1.0/self.interpol.n[self.k]) - self.r[self.k-1]**(3.0+1.0/self.interpol.n[self.k]))
            avg_shear += (-0.5*self.printparams.G*1.0/self.interpol.K[self.k])**(1.0/self.interpol.n[self.k])*1.0/(2.0+1.0/self.interpol.n[self.k])*(self.printparams.R**(2.0+1.0/self.interpol.n[self.k]) - self.r[self.k-1]**(2.0+1.0/self.interpol.n[self.k]))
            avg_visc  += self.interpol.K[self.k] * (-0.5*self.printparams.G*1.0/self.interpol.K[self.k])**(1.0-1.0/self.interpol.n[self.k])*1.0/(3.0-1.0/self.interpol.n[self.k])*(self.printparams.R**(3.0-1.0/self.interpol.n[self.k]) - self.r[self.k-1]**(3.0-1.0/self.interpol.n[self.k]))
            # Add the intermediate intervals
            for j in range(0,self.k,1) :
                flowrate += 0.5*self.c[j]*(self.r[j]**2 - self.r[j-1]**2 ) - (-0.5*self.printparams.G*1.0/self.interpol.K[j])**(1.0/self.interpol.n[j])*self.interpol.n[j]*1.0/(1.0+self.interpol.n[j])*1.0/(3.0+1.0/self.interpol.n[j])*(self.r[j]**(3.0+1.0/self.interpol.n[j]) - self.r[j-1]**(3.0+1.0/self.interpol.n[j]))
                avg_shear += (-0.5*self.printparams.G*1.0/self.interpol.K[j])**(1.0/self.interpol.n[j])*1.0/(2.0+1.0/self.interpol.n[j])*(self.r[j]**(2.0+1.0/self.interpol.n[j]) - self.r[j-1]**(2.0+1.0/self.interpol.n[j]))
                avg_visc  += self.interpol.K[j] * (-0.5*self.printparams.G*1.0/self.interpol.K[j])**(1.0-1.0/self.interpol.n[j])*1.0/(3.0-1.0/self.interpol.n[j])*(self.r[j]**(3.0-1.0/self.interpol.n[j]) - self.r[j-1]**(3.0-1.0/self.interpol.n[j]))                
        flowrate  *= 2.0 * math.pi
        avg_shear *= 2.0 * math.pi
        avg_visc  *= 2.0 * math.pi
        # Store in class members
        self.flowrate = flowrate
        self.avg_shearrate = avg_shear* 1.0/area
        self.avg_viscosity = avg_visc* 1.0/area
        self.avg_velocity = self.flowrate * 1.0/area           
    

    def calc_discrete_averages(self) :
        """Calculate the average velocity, shear rate and viscosity from the sampled data"""
        area      = math.pi * self.printparams.R * self.printparams.R
        flowrate  = 0
        avg_shear = 0
        avg_visc  = 0
        avg_vel   = 0
        for i in range(1,self.samples+1,1) :
            ringarea   = math.pi * (self.data[0][i]*self.data[0][i] - self.data[0][i-1]*self.data[0][i-1])
            flowrate  += ringarea * self.data[1][i-1]
            avg_shear += ringarea * self.data[2][i-1]
            avg_visc  += ringarea * self.data[3][i-1]
        avg_vel    = flowrate * 1.0/area
        avg_shear /= area
        avg_visc  /= area
        return avg_vel,avg_shear,avg_visc,flowrate        
        
    def calc_profile_data(self) : 
        """Calculate and store the velocity, shear rate and viscosity profiles"""
        # Determine the range in r (radial position)
        for i in self.samplerange :
            # Find corresponding interval to radial position
            r = i*1.0/self.samples * self.printparams.R
            interval = self.find_interval(r)
            # Radial position
            self.data[0][i] = r
            # Velocity
            self.data[1][i] = self.calc_velocity(r, interval)
            # Shear rate
            self.data[2][i] = self.calc_shearrate(r, interval)
            # Viscosity
            self.data[3][i] = self.interpol.calc_powerlaw_visc(self.calc_shearrate(r, interval),interval)
            # interpolation interval
            self.data[4][i] = interval
            

    def save_profiles(self, filename) : 
        """Save the calculated velocity, shear rate and viscosity profiles to a file"""
        datafile = open(filename, 'w')
        datafile.write("#r\tvelocity(m/s)\tshearrate(1/s)\tviscosity(Pa*s)\n")
        for r in self.samplerange :
            datafile.write("%.6E\t%.6E\t%.6E\t%.6E\t%i\n" % (self.data[0][r],self.data[1][r],self.data[2][r],self.data[3][r],self.data[4][r],))
            
            
    def save_averages(self, filename) :
        """Save the averaged quantities"""
        datafile = open(filename, 'w')
        datafile.write("Averaged quantities:\n")
        datafile.write("flowrate \t = %.6E m^3/s\n" % (self.flowrate))
        datafile.write("avg_vel \t = %.6E m/s\n" % (self.avg_velocity))
        datafile.write("avg_shear \t = %.6E 1/s\n" % (self.avg_shearrate))
        datafile.write("avg_visc \t = %.6E Pa*s\n" % (self.avg_viscosity))       


    def plot_velocity(self, draw_intervals=None, save_figure=None, ymin=None, ymax=None, xlabel=None, ylabel=None, title=None) : 
        """Plot the calculated velocity profiles and/or save the data to a file"""
        # Create plot
        fig = plt.figure(1, figsize=(12,9))
        ax = fig.add_subplot(111)
        # Add plot title
        if title == None :
            fig.suptitle('Velocity profile', fontsize=20)
        else :
            fig.suptitle(title, fontsize=20)
        # Add axis labels
        if xlabel == None :
            ax.set_xlabel(r'radial position $r$ in $\mathrm{m}$', fontsize=18)
        else : 
            ax.set_xlabel(xlabel, fontsize=18)
        if ylabel == None :
            ax.set_ylabel(r'velocity $u$ in $\frac{\mathrm{m}}{\mathrm{s}}$', fontsize=18) 
        else :
            ax.set_ylabel(ylabel, fontsize=18)
        # Set plot ranges
        ax.set_xlim(-self.printparams.R, self.printparams.R);
        if ymin == None and ymax == None :
            plotymin = 0.0
            plotymax = 1.1 * self.c[self.k]
        else :
            plotymin = ymin
            plotymax = ymax
        ax.set_ylim(plotymin, plotymax);
        # Draw the interpolated intervals
        if draw_intervals != None :
            # Dummy plot for interpolation intervals legend
            plt.plot([], linewidth=0.5, color='0.25', label=r'intervals: $r_i$')            
            for i in self.krange :
                ax.axvline(self.r[i], linewidth=0.5, color='0.25')
                ax.axvline(-self.r[i], linewidth=0.5, color='0.25')
        # Print the fluid profile
        ax.plot(self.data[0][self.samplerange], self.data[1][self.samplerange], 'b-', label='Piecewise solution')
        ax.plot(-self.data[0][self.samplerange], self.data[1][self.samplerange], 'b-') # Mirror image
        # Further plot options
        # Key options
        ax.legend()        
        # Invert x axis labels (only positive radial position, but mirrored for better visibility)
        ax.set_xticklabels([str(abs(x)) for x in ax.get_xticks()])
        # Show plot in Console    
        plt.show()        
        # Save the plot as PNG file
        if save_figure != None :
            fig.savefig(save_figure)
            
                
    def plot_shearrate(self, draw_intervals=None, save_figure=None, ymin=None, ymax=None, xlabel=None, ylabel=None, title=None) : 
        """Plot the calculated shear rate profiles"""
        # Create plot
        fig = plt.figure(1, figsize=(12,9))
        ax = fig.add_subplot(111)
        # Add plot title
        if title == None :
            fig.suptitle('Shear rate profile', fontsize=20)
        else :
            fig.suptitle(title, fontsize=20)
        # Add axis labels
        if xlabel == None :
            ax.set_xlabel(r'radial position $r$ in $\mathrm{m}$', fontsize=18)
        else : 
            ax.set_xlabel(xlabel, fontsize=18)
        if ylabel == None :
            ax.set_ylabel(r'shear rate $\dot{\gamma}$ in $\frac{1}{\mathrm{s}}$', fontsize=18)
        else :
            ax.set_ylabel(ylabel, fontsize=18)
        # Set plot ranges
        ax.set_xlim(-self.printparams.R, self.printparams.R);
        if ymin == None and ymax == None :
            plotymin = 0.0
            plotymax = 1.1 * self.calc_shearrate(self.printparams.R,self.k)
        else :
            plotymin = ymin
            plotymax = ymax
        ax.set_ylim(plotymin, plotymax);
        # Draw the interpolated intervals
        if draw_intervals != None :
            # Dummy plot for interpolation intervals legend
            plt.plot([], linewidth=0.5, color='0.25', label=r'intervals: $r_i$')            
            for i in self.krange :
                ax.axvline(self.r[i], linewidth=0.5, color='0.25')
                ax.axvline(-self.r[i], linewidth=0.5, color='0.25')
        # Print the fluid profile
        ax.plot(self.data[0][self.samplerange], self.data[2][self.samplerange], 'b-', label='Piecewise solution')
        ax.plot(-self.data[0][self.samplerange], self.data[2][self.samplerange], 'b-') # Mirror image
        # Further plot options
        # Key options
        ax.legend()        
        # Invert x axis labels (only positive radial position, but mirrored for better visibility)
        ax.set_xticklabels([str(abs(x)) for x in ax.get_xticks()])
        # Show plot in Console    
        plt.show()        
        # Save the plot as PNG file
        if save_figure != None :
            fig.savefig(save_figure)
            
            
    def plot_viscosity(self, draw_intervals=None, save_figure=None, ymin=None, ymax=None, xlabel=None, ylabel=None, title=None, logarithmic=None) : 
        """Plot the calculated viscosity profiles"""       
        # Create plot
        fig = plt.figure(1, figsize=(12,9))
        ax = fig.add_subplot(111)
        # Add plot title
        if title == None :
            fig.suptitle('Viscosity profile', fontsize=20)
        else :
            fig.suptitle(title, fontsize=20)
        # Add axis labels
        if xlabel == None :
            ax.set_xlabel(r'radial position $r$ in $\mathrm{m}$', fontsize=18)
        else : 
            ax.set_xlabel(xlabel, fontsize=18)
        if ylabel == None :
            ax.set_ylabel(r'dynamic viscosity $\eta$ in $\mathrm{Pa\,s}$', fontsize=18)
        else :
            ax.set_ylabel(ylabel, fontsize=18)
        # Set plot ranges
        ax.set_xlim(-self.printparams.R, self.printparams.R);
        if ymin == None and ymax == None :
            plotymin = 0.9 * self.calc_viscosity(self.printparams.R,self.k)
            plotymax = 1.1 * self.interpol.K[0]
        else :
            plotymin = ymin
            plotymax = ymax
        ax.set_ylim(plotymin, plotymax);
        # Draw the interpolated intervals
        if draw_intervals != None :
            # Dummy plot for interpolation intervals legend
            plt.plot([], linewidth=0.5, color='0.25', label=r'intervals: $r_i$')            
            for i in self.krange :
                ax.axvline(self.r[i], linewidth=0.5, color='0.25')
                ax.axvline(-self.r[i], linewidth=0.5, color='0.25')
        # Print the fluid profile
        if logarithmic != None :
            ax.semilogy(self.data[0][self.samplerange], self.data[3][self.samplerange], 'b-', label='Piecewise solution')
            ax.semilogy(-self.data[0][self.samplerange], self.data[3][self.samplerange], 'b-') # Mirror image
        else :
            ax.plot(self.data[0][self.samplerange], self.data[3][self.samplerange], 'b-', label='Piecewise solution')
            ax.plot(-self.data[0][self.samplerange], self.data[3][self.samplerange], 'b-') # Mirror image
        # Further plot options
        # Key options
        ax.legend()        
        # Invert x axis labels (only positive radial position, but mirrored for better visibility)
        ax.set_xticklabels([str(abs(x)) for x in ax.get_xticks()])
        # Show plot in Console    
        plt.show()        
        # Save the plot as PNG file
        if save_figure != None :
            fig.savefig(save_figure)