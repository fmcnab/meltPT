from meltPT import *
import yaml
import sys
import matplotlib.pyplot as plt

# ---- Read parameter file
if len(sys.argv) < 2:
    print("No parameter file provided, exiting.")
    print()
    print("The correct syntax is:")
    print("$ python3 run.py /path/to/paramter/file")
    sys.exit()
param_filename = sys.argv[1]
with open(param_filename, 'r') as f:
    params = yaml.safe_load(f)
    
# ---- Read and process data
s = Suite(
    params['data_file'], 
    src_FeIII_totFe=params['src_FeIII_totFe'], 
    min_MgO=params['min_MgO'], 
    min_SiO2=params['min_SiO2'])
s.backtrack_compositions()
s.compute_pressure_temperature()

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
max_P = s.PT['P'].max() + 1.
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# ---- Fit Tp to suite
if params['suite_Tp']['find']:
    s.find_suite_potential_temperature(mantle, find_bounds=True)
    if params['suite_Tp']['plot']:
        plt.plot(T_sol, P_sol, "k")
        lab=r"$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
            s.potential_temperature, 
            s.upper_potential_temperature - s.potential_temperature,
            s.potential_temperature - s.lower_potential_temperature)
        plt.plot(s.path.T, s.path.P, "--", label=lab)
        plt.plot(s.upper_path.T, s.upper_path.P, ":")
        plt.plot(s.lower_path.T, s.lower_path.P, ":")
        plt.scatter(s.PT['T'], s.PT['P'], marker="*")
        plt.xlabel(r"Temperature [$^\circ$C]")
        plt.ylabel("Pressure [GPa]")
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()