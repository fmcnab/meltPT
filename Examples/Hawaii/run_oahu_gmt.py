
import pandas as pd
import pygmt
import scipy
from meltPT import *
import matplotlib.pyplot as plt

# ----
# Generate Results
# ----

# ---- Oahu data
s = Suite("Oahu.csv", src_FeIII_totFe=0.15, min_MgO=8.5, min_SiO2=40.)
s.backtrack_compositions()
s.compute_pressure_temperature()


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
# lz.parameters['Mcpx'] = 0.15 # shorttle
# lz.CP = 1187. # shorttle
# lz.alphas = 30. # shorttle
# lz.DeltaS = 407. # shorttle
mantle = m.mantle([lz], [1], ['Lz'])
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]


# ---- Fit Tp to suite
s.find_suite_potential_temperature(mantle, find_bounds=True)


# ----
# Plot Results
# ----

# Load sample penguins data and convert 'species' column to categorical dtype
#df = pd.read_csv("https://github.com/mwaskom/seaborn-data/raw/master/penguins.csv")
#df.species = df.species.astype(dtype="category")

# Use pygmt.info to get region bounds (xmin, xmax, ymin, ymax)
# The below example will return a numpy array like [30.0, 60.0, 12.0, 22.0]
region = [1200., 1800., 0., 8.]

# Make a 2D categorical scatter plot, coloring each of the 3 species
# differently
fig = pygmt.Figure()

# Generate a basemap of 10 cm x 10 cm size
fig.basemap(
    region=region,
    projection="X10c/10c",
    frame=[
        'xafg+l"Temperature (oC)"',
        'yafg+l"Pressure (GPa)"',
    ],
)

# Define a colormap to be used for three categories, define the range of the
# new discrete CPT using series=(lowest_value, highest_value, interval),
# use color_model="+cAdelie,Chinstrap,Gentoo" to write the discrete color
# palette "inferno" in categorical format and add the species names as
# annotations for the colorbar
#pygmt.makecpt(cmap="inferno", series=(0, 2, 1), color_model="+cAdelie,Chinstrap,Gentoo")

fig.plot(
    # Use bill length and bill depth as x and y data input, respectively
    x=s.PT['T'],
    y=s.PT['P'],
    # Vary each symbol size according to another feature (body mass,
    # scaled by 7.5*10e-5)
    size=0.1,
    # Points colored by categorical number code
    color='k',
    # Use colormap created by makecpt
    #cmap=True,
    # Do not clip symbols that fall close to the map bounds
    no_clip=True,
    # Use circles as symbols with size in centimeter units
    style="cc",
    # Set transparency level for all symbols to deal with overplotting
    #transparency=40,
)

# Add colorbar legend
#fig.colorbar()

fig.show()
