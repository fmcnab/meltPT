# ---- Import libraries
from meltPT import *

# In this tutorial, we discuss how to implement your own reverse fractionation
# scheme in meltPT. As with the thermobarometers in Tutorial 6, this
# backtracking is carried out by a separate python class. For now, we only
# provide a scheme that corrects for the effects of olivine, but we hope to add
# more functionality in the future (contributions are welcome!).

# The key features of a backtracking class are:
# - contains a method called "backtrack_sample_composition"
# - "backtrack_sample_composition" should accept only a single-row pandas
#   dataframe containing the compositions in oxide wt% to be backtracked
# - "backtrack_sample_composition" should return a dictionary or series
#   containing the backtracked compositions. We expect these compositions to be
#   in oxide wt percent, and column headings should have the suffix
#   "_primary_wt"
# - any parameters or setting up needed prior to backtracking should be done
#   during instanciation. The class instance will then be passed to a Suite
#   instance for use.

# Let's design a simple minumum working example.
# This "backtracker" will simply add silica until some threshold is reached.
class ExampleBacktracker:
    
    def __init__(self, target_SiO2, dm):
        self.target_SiO2 = target_SiO2
        self.dm = dm
        
    def add_silica(self):
        
        # Add the specified amount of silica
        self.oxide_wt_hydrous['SiO2'] += self.dm
        
        # Renormalise to a total of 100% using meltPT's normalise function
        self.oxide_wt_hydrous = normalise(self.oxide_wt_hydrous)
        
    def backtrack_sample_composition(self, df):
        
        # First we will extract relevant major oxides from the input dataframe.
        # meltPT includes a list, MAJOR_OXIDES, containing each of the oxides
        # used by its thermobarometers. We will extract each of these
        # concentrations and assign them to oxide_wt_hydrous.
        self.oxide_wt_hydrous = {}
        for ox in MAJOR_OXIDES:
            self.oxide_wt_hydrous[ox] = df[ox]
        self.oxide_wt_hydrous = normalise(self.oxide_wt_hydrous)
        
        # Next let's iteratively add silica until we reach the specified
        # threshold.
        while self.oxide_wt_hydrous['SiO2'] < self.target_SiO2:
            self.add_silica()
            
        # Now we need to add some suffixes and return the new compositions.
        primary_oxide = {}
        for phase in self.oxide_wt_hydrous:
            primary_oxide[phase + "_primary_wt"] = self.oxide_wt_hydrous[phase]
            
        # And return.
        return primary_oxide


# Now let's try it out on our sample from Plank & Forsyth (2016)
s = Suite("../Data/PF16_UT09DV04.csv", src_FeIII_totFe=0.17, src_Fo=0.9)

# We will set up our backtracking class, with a target SiO2 content of 60% and
# an incremement of 0.1 wt%. Then we pass it to the Suite.backtrack_compositions
# method.
b = ExampleBacktracker(target_SiO2=60., dm=0.001)
s.backtrack_compositions(backtracker=b)

# Take a look at the result. In particular, we see that the SiO2 concentration
# is close to 60 wt%. All other phases now have reduced concentration.
print(s.primary)
print(s.primary['SiO2_primary_wt'])

# We can have a go at calculating equilibration presures and temperatures.
# The results will be nonsense of course.
s.compute_pressure(method="PF16")
print(s.PT)
