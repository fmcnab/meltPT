# ---- Importing
# Let's start by importing some libraries. If you have installed meltPT
# correctly, this should work!
from meltPT import *
import pyMelt as m
import matplotlib.pyplot as plt


# ---- Reading data and initialising the Suite object
# Start by importing data from your csv. In this case our data are stored in
# a file called "PF16_UT09DV04.csv". To keep this example simple, the input infile
# contains a single sample, UT09DV04, from Plank & Forsyth (2016). If necessary, 
# please edit the path within Suite so that the csv file can be found.
s = Suite("./Data/PF16_UT09DV04.csv", src_FeIII_totFe=0.17)

# You have now created an instance of meltPT's Suite class, containing the
# sample data. Data are stored in a pandas dataframe, which you can preview by
# typing:
print(s.data)


# ---- Backtrack compositions
# The next step is backtracking the sample's composition. The aim here is to 
# account for the effects of fractional crystallisation of olivine and estimate
# the sample's "primary" composition: i.e., its composition when it last
# equilibrated with the mantle. To do so, use the Suite's backtrack_compositions
# method. Here we have set the verbose flag to True, so program will print
# updates at each interation.
s.backtrack_compositions(Kd=0.3, verbose=True)

# As you can see, the sample started with a Forsterite number of c. 0.85. The
# program then added olivine in equilibrium with the melt, until, after adding
# c. 14% olivine, it reached a Forsterite number of 90%, which we assumed for
# the mantle source.

# You have now created a new dataframe within the suite class containing the
# sample's primary composition:
print(s.primary)


# ---- Compute pressures & temperatures
# Now we can calculate pressures and temperatures at which the sample was last
# in equilibrium with the mantle. Use Suite's compute_pressure_temperature
# method. 
s.compute_pressure_temperature(method="PF16")

# The results are stored in a new dataframe:
print(s.PT)

# The calculated pressure of 2.07 GPa and temperature of 1370 oC are the same 
# as those from Plank & Forsyth (2016, their Table S8), which is good!


# ---- Fit a melting path
# Next we would like to link our estimate equilibration pressure and temperature
# to a model geotherm. In meltPT, we make use of the pyMelt package to compute
# adiabatic decompression melting paths. For this example we use pyMelt's
# implientation of Katz et al.'s (2003) lherzolite melting model to set up the
# mantle object.
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])

# We can now pass our mantle object to Suite's find_individual_potential_temperatures
# method. This will take a few seconds!
s.find_individual_potential_temperatures(mantle)

# This will create a new dataframe, called individual_potential_temperatures,
# which contains information about the result and a copy of the best-fitting
# pyMelt path object:
print(s.individual_potential_temperatures)


# ---- Variable Kd
# In the above example we used a fixed value for the partition coefficient
# that describes partitioning between olivine and melt, Kd. There is also an
# option to allow this value to vary as a function of melt Mg#, after Tamura
# et al. (2020, JoP). Let's repeat the steps above with this option.

# Read the csv, create a new Suite object.
s_varKd = Suite("./Data/PF16_UT09DV04.csv", src_FeIII_totFe=0.17)

# Backtrack. Variable Kd is the default option, so we can simply miss out that
# flag. Notice that in this case, we need to add more olivine, about 17%, to
# reach our target forsterite content.
s_varKd.backtrack_compositions(verbose=True)

# Compute pressure and temperature. The result is hotter and deeper than the
# fixed Kd case.
s_varKd.compute_pressure_temperature(method="PF16")
print(s_varKd.PT)

# Find the best fitting potential temperature.
s_varKd.find_individual_potential_temperatures(mantle)


# ---- Plot
# Now let's make a nice plot of our result!

# Initialise seperate plots for the fixed and variable Kd axes.
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(14,5))

# Set up pressure array.
P = np.arange(1., 4.5, 0.01)

# Plot the solidus.
axs[0].plot(lz.TSolidus(P), P, c="0.75", label="Solidus")
axs[1].plot(lz.TSolidus(P), P, c="0.75", label="Solidus")

# Plot an adiabat corresponding to the best-fitting potential temperature
adiabat_label = r"Adiabat: $T_p$ = %i $^{\circ}$C" % s.individual_potential_temperatures.iloc[0]['Tp']
axs[0].plot(
    mantle.adiabat(
        P, 
        s.individual_potential_temperatures.iloc[0]['Tp']), 
    P, ":", label=adiabat_label)
adiabat_label = r"Adiabat: $T_p$ = %i $^{\circ}$C" % s_varKd.individual_potential_temperatures.iloc[0]['Tp']
axs[1].plot(
    mantle.adiabat(
        P, 
        s_varKd.individual_potential_temperatures.iloc[0]['Tp']), 
    P, ":", label=adiabat_label)
    
# Plot melt path corresponding to best-fitting potential temperature
melt_label = r"Melting path: $T_p$ = %i $^{\circ}$C" % s.individual_potential_temperatures.iloc[0]['Tp']
axs[0].plot(
    s.individual_potential_temperatures.iloc[0].path.T, 
    s.individual_potential_temperatures.iloc[0].path.P,
    label=melt_label)
melt_label = r"Melting path: $T_p$ = %i $^{\circ}$C" % s_varKd.individual_potential_temperatures.iloc[0]['Tp']
axs[1].plot(
    s_varKd.individual_potential_temperatures.iloc[0].path.T, 
    s_varKd.individual_potential_temperatures.iloc[0].path.P,
    label=melt_label)
    
# Plot our sample!
sample_label = r"UT09DV04, $F$ = %.2f%%" % (s.individual_potential_temperatures.iloc[0]['F']*100.)
axs[0].plot(s.PT['T'], s.PT['P'], "o", label=sample_label)
sample_label = r"UT09DV04, $F$ = %.2f%%" % (s_varKd.individual_potential_temperatures.iloc[0]['F']*100.)
axs[1].plot(s_varKd.PT['T'], s_varKd.PT['P'], "o", label=sample_label)

# Do some formatting and reveal
for ax in axs:
    ax.set_xlabel(r"Temperature [$^{\circ}$C]")
    ax.set_xlim(1300., 1500.)
    ax.set_ylim(1., 3.5)
    ax.legend()
    ax.invert_yaxis()
axs[0].set_ylabel(r"Pressure [GPa]")
plt.show()


# ---- Testing different thermobarometric schemes
# meltPT includes several different thermobarometric schemes, which are simple
# to interchange. Here we will quickly show you how to do this. In Tutorial 4,
# we test out the different schemes more systematically, and in Tutorial 6, you
# can see how to implement your own scheme.

# To see a list of the available options, we can look at the
# compute_presssure_temperature method's doc string. Type the following into
# the prompt:
# s.compute_pressure_temperature?

s = Suite("./Data/PF16_UT09DV04.csv", src_FeIII_totFe=0.17, param_co2=True)
s.backtrack_compositions(Kd=0.3, verbose=True)

# First let's define a list of the schemes we would like to use. We will take
# a subset for now.
schemes =  ["PF16", "L09", "P08", "BK21", "SD20"]
schemesT = ["B93", "P07_2", "P07_4", "HA15"]

# Loop over the schemes saving the result in a separate list.
PT = {}
T = {}
for scheme in schemes:
    s.compute_pressure_temperature(method=scheme)
    PT[scheme] = s.PT.copy()
for scheme in schemesT:
    s.compute_temperature(method=scheme, P=PT['PF16']['P'])
    T[scheme] = s.PT.copy()
    
# Plot
plt.plot(lz.TSolidus(P), P, c="0.75", label="Solidus")
for scheme in schemes:
    plt.errorbar(PT[scheme]['T'], PT[scheme]['P'], xerr=PT[scheme]['T_err'], yerr=PT[scheme]['P_err'], label=scheme, marker="o")
for scheme in schemesT:
    plt.errorbar(T[scheme]['T'], T[scheme]['P'], xerr=T[scheme]['T_err'], label=scheme, marker="^")

plt.legend()
plt.xlim(1300., 1550.)
plt.ylim(1., 4.5)
plt.xlabel(r"Temperature [$^{\circ}$C]")
plt.ylabel(r"Pressure [GPa]")
plt.gca().invert_yaxis()
plt.show()
    