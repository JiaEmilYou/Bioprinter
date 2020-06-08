#!/usr/bin/python

# ============================================================================
#            IMPORT PACKAGES
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
import math


# ============================================================================
#            TODO-LIST
# ============================================================================
# - include other analytical models (e.g. reduced CY of Ben Fabry) ?
# - include CY fit of provided gamma(shear) data (must)


# ============================================================================
#            CLASS DEFINITIONS
# ============================================================================
class Interpolation : 
    """Interpolation class interpolates given analytical function with continuous, piecewise defined, power-law functions"""

    def __init__(self, gamma0=None, gammaN=None, Ninterpol=None, gammaStart=None, gammaEnd=None, analytical=None, samples=10000) :
        """Create a new Interpolation with given or default parameters"""
        # Check if correct number of arguments is given, otherwise print instructions
        if gamma0 == None or gammaN == None  or Ninterpol == None or analytical == None : 
            self.print_usage()
        else :
            # Bounds of range of the interpolation functions (only power-law intervals)
            self.gamma0     = gamma0
            self.gammaN     = gammaN
            # Number of interpolation intervals
            self.N          = Ninterpol
            # Bounds of full range of viscosity interpolation (including Newtonian limits)
            if gammaStart == None or gammaEnd == None :
                self.gammaStart = gamma0 * 1.0e-1
                self.gammaEnd   = gammaN * 1.0e1
            else :  
                self.gammaStart = gammaStart
                self.gammaEnd   = gammaEnd
            # Lists storing the parameters of the interpolation functions
            self.n          = np.ones(self.N+2)
            self.K          = np.zeros(self.N+2)
            self.gammai     = np.zeros(self.N+3)
            # Index range for full interpolation (including Newtonian limits)
            self.fullrange  = range(0,self.N+2,1)
            # The analytical form to be interpolated
            self.analytical = analytical
            # The number of samples used for data saving and plotting
            self.samples     = samples
            self.samplerange = range(0,self.samples+1,1)
            self.data        = np.zeros((4, self.samples+1), dtype=np.float64) # shear rate, viscosity(interpolated), viscosity(analytical), interval
    
        
    @staticmethod
    def print_usage() :
        """Method used to print input instructions."""
        print("\n --- Interpolation class usage --- \n")
        print("Interpolation() class requires at least 4 parameters:")
        print(" - gamma0 :\t\t\t lower shear rate limit for interpolation")
        print(" - gammaN :\t\t\t upper shear rate limit for interpolation")
        print(" - Ninterpol :\t\t\t number of interpolation intervals")
        print(" - analytical :\t\t\t analytical viscosity model as instance of Analytical_viscosity() class")
        print(" - (optional) gammaStart :\t lower shear rate limit for plotting")
        print(" - (optional) gammaEnd :\t upper shear rate limit for plotting")
        print(" - (optional) samples :\t\t number of points calculated for plotting and data saving\n")
        print("The following methods must be called:")
        print(" - Interpolation(<params>) :\t\t with all necessary parameters")
        print(" - calculate_interpolation() :\t\t performs the interpolation \n")
        print("The following methods are helpful:")
        print(" - save_interpolation(filename) : \t\t save the viscosity-shearrate data to a file")
        print(" - save_interpolation_parameters(filename) :\t save the interpolation parameters to a file")
        print(" - plot_interpolation(options) :\t\t plot the viscosity-shearrate relationship")
        print("   Options are:")
        print("   - plot_interpolation=<int> : \t\t plot the interpolated viscosity")
        print("   - plot_analytical=<int> : \t\t\t plot the analytical viscosity (numbers give the order in the plot)")
        print("   - draw_intervals=<int> : \t\t\t indicate the interval boundaries with lines")
        print("   - draw_last_interval=<int> : \t\t indicate the last interval used for profile calculation with a thick line (provide last interval with Profiles().get_last_interval())")
        print("   - save_figure=<filename> : \t\t\t save the plot to a file\n\n")

    
    
    def set_interpolation_range(self, gamma0, gammaN) :
        """Set the interpolation range bounds"""
        self.gamma0 = gamma0
        self.gammaN = gammaN
        
        
    def set_full_range(self, gammaStart, gammaEnd) :
        """Set the bounds of the full gamma range, which includes the interpolation range"""
        self.gammaStart = gammaStart
        self.gammaEnd   = gammaEnd
    
    
    def set_interpolation_intervals(self, Nintervals) :
        """Set the number of interpolation intervals in the interpolation range"""
        self.N = Nintervals
        self.reinit_params()
        
        
    def set_analytical_viscostiy(self, analytical) :
        """Set the analytical viscosity model"""
        self.analytical = analytical
        
        
    def reinit_params(self) :
        """Re-initializes the parameter arrays, if number of interpolation intervals, N, changes"""
        self.n          = np.ones(self.N+2)
        self.K          = np.zeros(self.N+2)
        self.gammai     = np.zeros(self.N+3)
        self.fullrange  = range(0,self.N+2,1)        


    def calc_gammai(self, i) :
        """Calculates the upper bound of interval, assuming log-equidistant partitioning of the range"""
        # Equidistant in log scale
        gammai = self.gamma0 * (self.gammaN *1.0/self.gamma0)**(1.0*i * 1.0/self.N)
        return gammai


    def calc_ni(self, i) :
        """Calculates the power-law exponent of interpolation interval i"""
        # The 0th and (n+1)th interval are Newtonian, the exponent equals 1
        if i==0 : 
            return 1.0
        elif i==self.N+1 :
            return 1.0
        else:
            # The exponents of the other intervals are calculated using the range and the analytical form
            ni = 1.0 - self.N * 1.0/(math.log(self.gammaN * 1.0/self.gamma0)) * math.log( self.analytical.calc_visc(self.gammai[i-1]) / self.analytical.calc_visc(self.gammai[i]) )
            return ni
 
   
    def calc_Ki(self,i) :
        """Calculates the power-law consistency parameter of interpolation interval <i>"""
        # The Newtonian regions have a constant viscosity, corresponding to the viscosity
        # at the beginning of the first interval (i=1) and
        # at the end of the last interval (i=N)
        if i==0 :
            return self.analytical.calc_visc(self.gamma0)
        elif i==self.N+1 :
            return self.analytical.calc_visc(self.gammaN)
        else :
            Ki = math.sqrt(self.analytical.calc_visc(self.gammai[i-1]) * self.analytical.calc_visc(self.gammai[i])) * (self.gammai[i-1] * self.gammai[i])**(0.5*(1.0-self.n[i]))
            return Ki
          
          
    def calc_powerlaw_visc(self, gamma, i) :
        """Calculate the interpolated viscosity at shear rate <gamma> in interval <i>"""
        visc = self.K[i] * gamma**(self.n[i] - 1.0)
        return visc

    def find_interval(self,gamma) :
        """Determine the index of the interval that includes the shear rate <gamma>"""
        i = 0
        while( self.gammai[i] < gamma) :
            i += 1
        return i
   
    def calculate_interpolation(self) :
        """Fills the parameter lists for the interpolation function parameters"""
        # Each list is filled in a separate for-loop, because the calculation requires the previously calculated values       
        # First, the complete shear rate (gamma) range is defined, which includes the interpolated range and the Newtonian regions
        self.gammai[-1]       = self.gammaStart # gammai[-1] equals gammai[N+2]
        self.gammai[self.N+1] = self.gammaEnd
        for i in range(0,self.N+1,1) :
            self.gammai[i]    = self.calc_gammai(i)         
        # Calculate the power-law exponents and consistency parameters of the interpolation functions
        # This range includes the Newtonian regions
        for i in self.fullrange :
            self.n[i] = self.calc_ni(i)
            
        for i in self.fullrange :
            self.K[i] = self.calc_Ki(i)
        # Store data in array
        self.store_data()

        
    def store_data(self) : 
        """Store the plotting data for viscosity(shear rate) in data array"""
        # Determine the range in r (radial position)
        for i in self.samplerange :
            # Find corresponding interval to shear rate gamma (equidistant on logarithmic scale)
            gamma = self.gammaStart * (self.gammaEnd*1.0/self.gammaStart)**(i*1.0/self.samples)
            interval = int(self.find_interval(gamma))
            # Shear rate
            self.data[0][i] = gamma
            # Viscosity (interpolated)
            self.data[1][i] = self.calc_powerlaw_visc(gamma, interval)
            # Viscosity (analytical)
            self.data[2][i] = self.analytical.calc_visc(gamma)
            # interpolation interval
            self.data[3][i] = interval
        
        
    def save_interpolation(self, filename) :
        """Save the interpolated and analytical viscosity as function of the shear rate"""
        if filename == None :
            print('Provide filename and path: save_interpolation("path/to/data.dat")')
        else :
            filehandle = open(filename, 'w')
            filehandle.write("# shearrate\tvisc(interpol)\tvisc(analytical)\tinterval")
            for i in self.samplerange :
                filehandle.write("%.6E\t%.6E\t%.6E\t%i\n" % (self.data[0][i],self.data[1][i],self.data[2][i],self.data[3][i]) )


    def save_interpolation_parameters(self, filename) :
        """Save the interpolation function parameters of each interval to a file"""
        if filename == None :
            print('Provide filename and path: save_interpolation_params("path/to/data.dat")')
        else :        
            filehandle = open(filename, 'w')
            filehandle.write("# i\tgammai\t\tKi\t\tni\n")
            for i in self.fullrange :
                filehandle.write("%i\t%.6E\t%.6E\t%.6E\n" % (i, self.gammai[i], self.K[i], self.n[i]))

    
    def plot_interpolation(self, plot_interpolation=1, plot_analytical=None, draw_intervals=None, draw_last_interval=None, save_figure=None, ymin=None, ymax=None, xlabel=None, ylabel=None, title=None) :
        """Plot the interpolation function (and optional the analytical form) in a log-log plot"""
        # Create plot
        fig = plt.figure(1, figsize=(12,9))
        fig.add_subplot(111)
        # Add plot title
        if title == None :
            fig.suptitle('Viscosity interpolation', fontsize=20)
        else :
            fig.suptitle(title, fontsize=20)
        # Add axis labels
        if xlabel == None :
            plt.xlabel(r'shear rate $\dot{\gamma}$ in $\frac{1}{\mathrm{s}}$', fontsize=18)
        else : 
            plt.xlabel(xlabel, fontsize=18)
        if ylabel == None :
            plt.ylabel(r'Dynamic viscosity $\eta(\dot{\gamma})$ in $\mathrm{Pa\,s}$', fontsize=18)
        else :
            plt.ylabel(ylabel, fontsize=18)

        # Set plot ranges
        plt.xlim(self.gammaStart, self.gammaEnd);
        if ymin != None and ymax != None : 
            plot_ymin = ymin
            plot_ymax = ymax
        else :
            plot_ymin = 1.0e-1*self.K[self.N+1]
            plot_ymax = 1.0e1*self.K[0]           
        plt.ylim(plot_ymin, plot_ymax);  

        # Draw the interpolated intervals
        if draw_intervals != None :
            # Dummy plot for interpolation intervals legend
            plt.loglog([], linewidth=0.5, color='0.25', label=r'Interpolation intervals: $\dot{\Gamma}_i$')
            # Draw intervals as vertical lines
            for i in self.fullrange :
                plt.axvline(self.gammai[i], linewidth=0.5, color='0.25')
        
        # Draw the k-th interval from flow profile calculation
        if draw_last_interval != None : 
            plt.axvline(self.gammai[draw_last_interval], linewidth=2.0, color='green', label=r'Last interval of profile calculation: $\dot{\Gamma}_{k=%i} = %.3E}$'%(draw_last_interval, self.gammai[draw_last_interval]))
        
        # Determine the plot order: interpolation in front of or behind analytical
        if plot_analytical != None :
            if plot_interpolation < plot_analytical :
                # Print the analytical viscosity model
                plt.loglog(self.data[0][self.samplerange], self.data[2][self.samplerange], 'r-', label=r'Carreau-Yasuda model: $\tilde{\eta}(\dot{\gamma})$') # : $\tilde{\eta}(\dot{\gamma}) = \eta_\infty + \frac{\eta_o - \eta_\infty}{(1 + (K\dot{\gamma})^{a_1})^{a_2(a_1)}}$            
                # Print the interpolation function
                plt.loglog(self.data[0][self.samplerange], self.data[1][self.samplerange], 'b-', label=r'Interpolation: $\eta(\dot{\gamma})$') # : $\eta(\dot{\gamma})$       
            else :
                # Print the interpolation function
                plt.loglog(self.data[0][self.samplerange], self.data[1][self.samplerange], 'b-', label=r'Interpolation: $\eta(\dot{\gamma})$') # : $\eta(\dot{\gamma})$       
                # Print the analytical viscosity model
                plt.loglog(self.data[0][self.samplerange], self.data[2][self.samplerange], 'r-', label=r'Carreau-Yasuda model: $\tilde{\eta}(\dot{\gamma})$') # : $\tilde{\eta}(\dot{\gamma}) = \eta_\infty + \frac{\eta_o - \eta_\infty}{(1 + (K\dot{\gamma})^{a_1})^{a_2(a_1)}}$
        else : 
            # Print only the interpolation function
            plt.loglog(self.data[0][self.samplerange], self.data[1][self.samplerange], 'b-', label=r'Interpolation: $\eta(\dot{\gamma})$') # : $\eta(\dot{\gamma})$                   

        # Key options
        plt.legend(loc='lower left')
        
        # Show plot in Console    
        plt.show()
        
        # Print plot to file
        if save_figure != None :
            fig.savefig(save_figure)            
        
  

class Analytical_Viscosity :    
    """Analytical_viscosity class holds the parameters for the analytical viscosity model"""
    def __init__(self, eta0=None, etainf=None, K=None, a1=None, a2=None) :
        """Initialize the parameters of the Carreau-Yasuda model"""
        # Check for correct number of arguments
        if eta0 == None or etainf == None or K == None or a1 == None or a2 == None :
            self.print_usage()
        else :
            self.eta0   = eta0
            self.etainf = etainf
            self.K      = K
            self.a1     = a1
            self.a2     = a2
   
    @staticmethod     
    def print_usage() :
        """Method used to print input instructions."""
        print(" --- Analytical_Viscosity class usage ---")
        print("\nAnalytical_Viscosity() class requires 5 parameters (Carreau-Yasuda viscosity model):")
        print(" - eta0 :\t viscosity in the limit of zero shear rate")
        print(" - etainf :\t viscosity in the limit of infinite shear rate")
        print(" - K :\t\t consistency index or inverse of corner shear rate")
        print(" - a1 :\t\t first exponent")
        print(" - a2 :\t\t second exponent\n\n")
        
    
    def set_parameters(self, eta0, etainf, K, a1, a2) :
        """Set the parameters of the Carreau-Yasuda model"""
        self.eta0   = eta0
        self.etainf = etainf
        self.K      = K
        self.a1     = a1
        self.a2     = a2
        
        
    def calc_visc(self, gamma) :
        """Calculate the viscosity at given shear rate for the analytical viscosity model"""
        visc = self.etainf + (self.eta0 - self.etainf) / ( (1.0 + (self.K * gamma)**(self.a1) )**(self.a2*1.0/self.a1) )
        return visc
        
        
    def save_parameters(self, filename, data_format='%.6E') :
        """Save the data of the analytical viscosity model to a file"""
        datafile = open(filename, 'w')
        datafile.write(("eta0\t = \t"+data_format+"\n") % self.eta0)
        datafile.write(("etainf\t = \t"+data_format+"\n") % self.etainf)
        datafile.write(("K\t = \t"+data_format+"\n") % self.K)
        datafile.write(("a1\t = \t"+data_format+"\n") % self.a1)
        datafile.write(("a2\t = \t"+data_format+"\n") % self.a2)
        datafile.write("\n# The viscosity is calculated from the shear rate (shear) using the formula: etainf + (eta0-etainf)/((1.0+(K*shear)^a1)^(a2/a1))\n")