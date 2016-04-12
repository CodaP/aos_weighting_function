# This python script handles the tasks of of loading and
# interpolating data from a sounding filed downloaded from the Wyoming
# site.  Your job will be to add the parts that do the actual
# radiative transfer.
# 
# Before making your own changes, you should take the time to
# carefully work through the code that is already here, study the
# syntax, and understand how it works.  Once you do, you will have a
# working introductory knowledge of python.


import matplotlib.pyplot as plt
from numpy import log, cos, exp, arange
import numpy as np
import gasabsr98
import pandas as pd

#=================

# Define function to compute saturation vapor pressure in Pa, given temperature in C.
def esat(TC):
  return(611.2*exp(17.67*TC/(TC+243.5)))

#=================

# Define function to compute saturation mixing ratio, given temperature in K and pressure in Pa
def wsat(TK,p):
  es = esat(TK-273.15)
  return(0.622*es/(p-es))

#=================
class Level:
  'Data for a single level in a sounding.  Inputs are pressure (Pa), temperature (C), dewpoint (C) '

  def __init__(self, Pres,  TempC,  DewPt):  # This internal method is
                                             # automatically executed
                                             # whenever a new object
                                             # of type 'Level' is
                                             # instantiated.
    self.Pres  = Pres                 #Pressure stored in Pascals
    self.TempC = TempC                #Temperature and dewpoint in C
    self.DewPt = DewPt
    self.TempK = self.TempC + 273.15  #Also save temperature in K
    self.Qvap = 0.622*esat(self.DewPt)/self.Pres  #Compute specific humidity

  def prt(self):
    print "%10.1f %8.1f %6.1f %6.1f %10.7f" % (self.Pres,self.TempC,self.DewPt,self.TempK,self.Qvap)


#=================
class Layer:
  'Layer of the atmosphere bounded by two successive levels in a sounding'

  def __init__(self,level1,level2):
    self.pbar = 0.5*(level1.Pres + level2.Pres)
    self.mass = (level1.Pres - level2.Pres)/9.80665
    self.qbar = 0.5*(level1.Qvap + level2.Qvap)
    self.Tbar = 0.5*(level1.TempK + level2.TempK)
    self.wvmass = self.qbar*self.mass
    self.DZ = 287.06*self.Tbar*(1.+0.61*self.qbar)*log(level1.Pres/level2.Pres)/9.80665
    self.rho = self.mass/self.DZ
    self.lwp = 0.0
    self.rhowv = self.wvmass/self.DZ

  def prt(self):
    print "%8.1f %6.1f %8.1f %6.1f %10.6f %8.2f %8.4f %10.6f " % (self.Zbot,self.DZ,self.pbar,self.Tbar,self.qbar,self.mass,self.wvmass,self.rho)


# specify radiosonde input file name
fname = "91408_2011_07_19_12.txt"

# Specify the number of leading lines to skip before treating lines as
# potentially valid data:
def read_sonde(fname, Nskip=7):

    # open the file for reading
    file = open(fname, "r")

    print "Successfully opened file ", file.name

    # read lines to be skipped 
    for i in range(Nskip):
      line = file.readline()

    # load the rest of file into a list of lines 
    lines = file.readlines()

    # close file
    file.close()


    # This object will store a list of all valid levels in the input sounding
    Sounding = []

    # We will use this variable to retain the value of the previous
    # pressure so as to ensure that we only keep lines for which the
    # pressure actually decreases relative to the previous line.

    pold = 100000.0

    # loop over input lines, extracting data elements that we care about
    # and storing them in our list of level data.

    for line in lines:
        words = line.split()   # From each line, create a list of 'words'
        if len(words) != 11: continue  # Skip if there is missing date in line
        p = float(words[0])    # First 'word' gives the pressure in hPa
        if (p < pold):         # Verify that pressure actually decreased relative to previous line
          t = float(words[2])  # Temperature in C
          dp =float(words[3])  # Dewpoint
          Sounding.append(Level(p*100, t, dp))  #Create level structure and append to list
          pold = p             # Save current pressure for future reference

    # Print the total number of levels read
    print len(Sounding), ' levels in raw input sounding'

    # Add a bogus level at 1 hPa to ensure that most of the mass in the
    # atmosphere is accounted for by the total profile.

    Sounding.append(Level(100.0, -50.0, -80.0))

    return Sounding

def interpolate_sounding(Sounding):

    # Now interpolate new levels into the sounding wherever gap in either
    # pressure or temperature is too large.  This is done in multiple
    # passes, the loop being terminated when no more gaps are detected.

    Sounding = list(Sounding)

    gaps = 1    # Ensure at least one pass
    while gaps:
      gaps = 0    # Now reset to no gaps unless we actually encounter one
      for i in range(len(Sounding)-1):
        if Sounding[i+1].Pres/Sounding[i].Pres < 0.9 or abs(Sounding[i+1].TempC - Sounding[i].TempC) > 1.0:
          gaps = 1
          pbar = 0.5*(Sounding[i].Pres + Sounding[i+1].Pres)
          tbar = 0.5*(Sounding[i].TempC + Sounding[i+1].TempC)
          dpbar = 0.5*(Sounding[i].DewPt + Sounding[i+1].DewPt)
          Sounding.insert(i+1,Level(pbar,tbar,dpbar))  #insert new interpolated level between the two that were compared

    print len(Sounding), ' levels in interpolated sounding'

    print '    p [Pa]    T [C] Td [C]  T [K]     q []'
     # List all of the levels in the raw sounding
    #for level in Sounding:
        #level.prt()

    return Sounding

def plot_sounding(Sounding):
    # Plot the interpolated sounding

    pres = []
    temp = []
    dwpt = []
    for level in Sounding:
      pres.append(level.Pres/100.0)
      temp.append(level.TempC)
      dwpt.append(level.DewPt)

    title = "Sounding from " + fname
    plt.title(title)
    plt.plot(temp, pres)
    plt.plot(dwpt,pres)
    plt.semilogy()
    plt.xlabel("Temperature [$^\circ$C]")
    plt.ylabel("Pressure [hPa]")
    plt.axis([-100.0, 40.0, 1000.0, 10.0])
    plt.savefig("Sounding.pdf")
    plt.savefig("Sounding.png")

    # clear figure and axes in prep for later plots
    plt.clf()
    plt.cla()

def create_atmosphere(Sounding):
    Atmosphere = []  # This object will contain a list of all layers in the sounding
    TPW = 0.0        # Total precipitable water vapor
    for i in range(len(Sounding)-1):
      Atmosphere.append(Layer(Sounding[i],Sounding[i+1]))
      TPW += Atmosphere[-1].wvmass

    print len(Atmosphere), ' layers in model atmosphere'
    print 'Total precipitable water = ',TPW,' kg/m^2'
    Atmosphere =  (pd.DataFrame(
        [{key:getattr(x, key) for key in ['pbar','Tbar','rho','rhowv','DZ','qbar','mass','wvmass']}
         for x in Atmosphere])
         .set_index('pbar', drop=False).sort_index())
    Atmosphere['Zbar'] = Atmosphere.set_index('pbar').sort_index(ascending=False).DZ.cumsum()
    return Atmosphere

# for layer in Atmosphere:
#   layer.prt()
# 
# print

# ================================================================ 
# WE ARE FINISHED LOADING AND PRE-PROCESSING OUR ATMOSPHERIC PROFILE.
# EVERYTHING BELOW THIS LINE PERTAINS TO THE FREQUENCY-DEPENDENT
# RADIATIVE TRANSFER CALCULATIONS
#=================================================================

def get_weighting_function_upwelling(f, layers, mu=1.):
    import gasabsr98
    layers['ke_air'] = layers.apply(lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[0], axis=1)
    layers['ke_wv'] = layers.apply(lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[1], axis=1)

    layers['tau'] = (layers.ke_air*(layers.mass - layers.wvmass) + layers.ke_wv*layers.wvmass)
    layers = layers.reset_index(drop=True)
    layers['tx_from_top'] = (layers
                    # Take tau(z)
                     .sort_values('Zbar', ascending=False).tau
                    # Compute total optical depths from toa to z for each z
                     .cumsum()
                    # Compute tx = e^(-tau/mu)
                    .apply(lambda tau: np.exp(-tau/mu)))

    # Compute the contribution from each layer
    layers['layer_weight'] = -(layers.sort_values('Zbar', ascending=False).tx_from_top.diff())

    # Get weight for each layer in m^-1
    layers['weighting_function'] = layers.layer_weight / layers.sort_values('Zbar', ascending=True).Zbar.diff() * 1000
    return layers

def get_weighting_function_downwelling(f, layers, mu=1):
    import gasabsr98
    layers['ke_air'] = layers.apply(lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[0], axis=1)
    layers['ke_wv'] = layers.apply(lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[1], axis=1)

    layers['tau'] = (layers.ke_air*(layers.mass - layers.wvmass) + layers.ke_wv*layers.wvmass)
    layers = layers.reset_index(drop=True)
    layers['tx_from_bottom'] = (layers
                    # Take tau(z)
                     .sort_values('Zbar', ascending=True).tau
                    # Compute total optical depths from toa to z for each z
                     .cumsum()
                    # Compute tx = e^(-tau/mu)
                    .apply(lambda tau: np.exp(-tau/mu)))

    # Compute the contribution from each layer
    layers['layer_weight'] = (layers.sort_values('Zbar', ascending=False).tx_from_bottom.diff())

    # Get weight for each layer in m^-1
    layers['weighting_function'] = layers.layer_weight / -layers.sort_values('Zbar', ascending=False).Zbar.diff() * 1000
    return layers

def get_weighting_function_satellite(f, layers, surf_emis, mu=1):
    # Compute reflectivity
    refl = 1.0-surf_emis

    down = get_weighting_function_downwelling(f, layers, mu).set_index('Zbar').ix[:, ('weighting_function','layer_weight','tx_from_bottom')]
    up = get_weighting_function_upwelling(f, layers, mu).set_index('Zbar').ix[:,('weighting_function', 'layer_weight','tx_from_top')]

    # Compose the two weighting functions
    both = down*refl*up.tx_from_top.min() + up

    # Include transmission in data structure for convenience
    both['tx_from_top'] = up.tx_from_top
    both['tx_from_bottom'] = down.tx_from_bottom
    return both


def bt(f, layers, stemp, surf_emis, mu=1.0):
    """
    Compute brightness temperature
    """
    # First get weighting function
    wf = get_weighting_function_satellite(f, layers, surf_emis, mu=mu)
    # Weighted-sum the temperatures and add surface emission
    return wf.tx_from_top.min()*stemp*surf_emis + sum(wf.layer_weight * layers.set_index('Zbar').Tbar)


### All plotting and usage of these functions was documented in an ipython notebook 
