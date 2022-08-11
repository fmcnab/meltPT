from meltPT import *

# In this tutorial we will show you how to apply your own thermobarometer
# in meltPT. Each thermobarometer is implemented as a separate python class.
# There are a few key features that these classes much contain to work properly.
# On top of these necessary features, you can include anything you need to get
# the thermobarometer working.

# The key features are:
# - accepts a single-row pandas dataframe on instantiation
# - includes one or more of the methods: compute_pressure_temperature,
#   compute_pressure, compute_temperature
# - these methods should return a dictionary containing calculated pressures
#   and temperatures and their uncertainties

# Let's design a simple minumum working example.
# No matter what you give it, it will return a pressure of 1 +/- 0.1 GPa and a
# temperature of 1000 +/- 100 oC. 
class ExampleThermobarometer:
    
    def __init__(self, df):
        self.df = df
        self.P_err = 0.1
        self.T_err = 100.
        
    def compute_pressure_temperature(self):
        return {'P': 1., 'P_err': self.P_err, 'T': 1000., 'T_err': self.T_err}
        
    def compute_pressure(self, T):
        return {'P': 1., 'P_err': self.P_err, 'T': T, 'T_err': np.nan}
        
    def compute_temperature(self, P):
        return {'P': P, 'P_err': np.nan, 'T': 1000., 'T_err': self.T_err}
        

# Now let's try it out on our sample from Plank & Forsyth (2016)
s = Suite("PF16_UT09DV04.csv", src_FeIII_totFe=0.17)
s.backtrack_compositions(target_Fo=0.9)

# Call the compute_pressure_temperature method and feed it our custom
# thermobarometer.
s.compute_pressure_temperature(method=ExampleThermobarometer)
print(s.PT)

# Now let's try out it out as a barometer.
s.compute_pressure(method=ExampleThermobarometer, T=1300.)
print(s.PT)

# ... and as a thermometer.
s.compute_temperature(method=ExampleThermobarometer, P=2.)
print(s.PT)