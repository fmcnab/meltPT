from meltPT import *
import matplotlib.pyplot as plt

##############
# Find PT of each ridge sample
##############


# ---- Gale et al., (2013) ridge data
s = Suite("Gale_Dataset_No_Backarc.csv", src_FeIII_totFe=0.12, min_SiO2=40., min_MgO=7.5)

s.backtrack_compositions()
s.compute_pressure_temperature()

# ---- Write CSV
s.write_to_csv("Ridge_output.csv")


##############
# Find Tp of each ridge segment
##############


# Input Ridge Database
df = parse_csv('Gale_Dataset_No_Backarc.csv', Ce_to_H2O=200., 
               src_FeIII_totFe=0.12, min_SiO2=40., min_MgO=7.5)

# Make array of all ridge segments
ridges = df['Province'].unique()

# ---- Set up output file
df_out = pd.DataFrame(columns = ['Ridge', 'Longitude', 'Latitude', 'N_Samples',
                                 'Tp', 'upper_Tp', 'lower_Tp'])

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 6., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]
path = mantle.adiabaticMelt(1330., Pstart=6., steps=101)

for i in range(0, (len(ridges) - 1), 1):
    
    # Make csv file for ridge segment
    ridge_df = df.loc[(df['Province'] == ridges[i])]
    ridge_df.to_csv('ridge.csv')

    # Parse ridge segment data to s and filter for key criteria
    s = Suite("ridge.csv", src_FeIII_totFe=0.12, min_SiO2=40., min_MgO=7.5)
    
    # Make sure ridge segment has data
    if len(s.data['Province']) > 0:

        # ---- Find PT of each sample
        s.backtrack_compositions()
        s.compute_pressure_temperature()  

        # ---- Fit Tp to suite
        s.find_suite_potential_temperature(mantle, find_bounds=True)
        plt.plot(T_sol, P_sol, "k")
        lab=r"$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
            s.potential_temperature, 
            s.upper_potential_temperature - s.potential_temperature,
            s.potential_temperature - s.lower_potential_temperature)

        # ---- Output ridge segment data to df_out
        df_out = df_out.append({'Ridge' : s.data['Province'][0],
                        'Longitude' : s.data['Longitude'].mean(),
                        'Latitude' : s.data['Latitude'].mean(),
                        'N_Samples' : len(s.data['Latitude']),
                        'Tp' : s.potential_temperature,
                        'upper_Tp' : s.upper_potential_temperature,
                        'lower_Tp' : s.lower_potential_temperature}, 
                        ignore_index = True)

# ---- Print df_out to output file
df_out.to_csv('ridges_Tp.csv')
