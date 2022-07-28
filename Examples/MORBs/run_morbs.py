from meltPT import *
import matplotlib.pyplot as plt

from meltPT import *
import matplotlib.pyplot as plt

#############
#Find PT of each ridge sample (MORB Figure a)
#############
print("Make MORB Figure (a):")
# ---- Read in Data
df = pd.read_csv("Gale_Dataset_No_Backarc.csv", sep=',')
df1 = df.loc[df['Latitude']>70.]
df = df.loc[df['Latitude']<60.]
df = df.loc[df['FeO']>0.]
df = pd.concat([df, df1], ignore_index=True)
df.to_csv("ridge_in.csv", sep=',')

print("Reading In Data...")
# ---- Gale et al., (2013) ridge data
s = Suite("ridge_in.csv", src_FeIII_totFe=0.14, Ce_to_H2O=200.,
          min_SiO2=40., min_MgO=8.5)

print("Backtracking Compositions...")
# ---- Backtrack Compositions
s.backtrack_compositions()

print("Calculating PT of Each Sample...")
# ---- Calculate PT
s.compute_pressure_temperature()


# ---- Plot PT Results For Each Sample
print("Plotting MORB Figure (a)...")
plt.scatter(s.PT['T'], s.PT['P'], marker="o", facecolors='none', 
            edgecolors='k', label="Sample")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.gca().invert_yaxis()
plt.title("MORB Figure (a)")

# # ---- Write CSV
# s.write_to_csv("Ridge_output.csv")
print("Done!")
print("---------------")
print()

#############
# Find Tp of each ridge segment
#############
print("Make MORB Figure (b):")

print("Read in Data...")
# Input Ridge Database
df = parse_csv('ridge_in.csv', Ce_to_H2O=200., 
                src_FeIII_totFe=0.14, min_SiO2=40., min_MgO=8.5)

print("Set Up Mantle Composition...")
# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 5., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

print("Loop through Ridges Calculating Tp...")
# Make array of all ridge segments
ridges = df['Province'].unique()

# ---- Set up output file
df_out = pd.DataFrame(columns = ['Ridge', 'Longitude', 'Latitude', 'N_Samples',
                                  'Tp', 'upper_Tp', 'lower_Tp'])

for i in range(0, (len(ridges) - 1), 1):
#for i in range(0, 1, 1):    
    # Make csv file for ridge segment
    ridge_df = df.loc[(df['Province'] == ridges[i])]
    ridge_df.to_csv('ridge.csv')

    # Parse ridge segment data to s and filter for key criteria
    s2 = Suite("ridge.csv", Ce_to_H2O=200., src_FeIII_totFe=0.14, 
              min_SiO2=40., min_MgO=8.5)
    
    # Make sure ridge segment has data
    if len(s.data['Province']) > 5:

        # ---- Find PT of each sample
        s2.backtrack_compositions()
        s2.compute_pressure_temperature()  

        if len(s2.PT['T']) > 5:

            # ---- Fit Tp to suite
            s2.find_suite_potential_temperature(mantle, find_bounds=True)
            print("Best-fitting melting model for Ridge",ridges[i],": Tp = %i oC." 
                  % (s2.potential_temperature))
            
            # ---- Output ridge segment data to df_out
            
            dict_result = {'Ridge' : [s2.data['Province'][0]],
                            'Longitude' : [s2.data['Longitude'].mean()],
                            'Latitude' : [s2.data['Latitude'].mean()],
                            'N_Samples' : [len(s2.data['Latitude'])],
                            'Tp' : [s2.potential_temperature],
                            'upper_Tp' : [s2.upper_potential_temperature],
                            'lower_Tp' : [s2.lower_potential_temperature]}
            df_result = pd.DataFrame.from_dict(dict_result)
            df_out = pd.concat([df_out, df_result], ignore_index = True)

print("Plot MORB Figure (b)")
# ---- Plot Histogram of Tp Results
bin_array=[1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575]
plt.hist(df_out['Tp'], bins=bin_array, histtype='bar', align='mid', orientation='vertical')
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Number of Ridges")
plt.title("MORB Figure (b)")
print("Done!")
print("---------------")
print()
# # ---- Print df_out to output file
# df_out.to_csv('ridges_Tp.csv')
