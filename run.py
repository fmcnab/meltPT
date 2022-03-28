from meltPT import *
import yaml
import sys
# import matplotlib
# matplotlib.use('TkAgg')
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
    src_FeIII_totFe=params['backtrack']['src_FeIII_totFe'], 
    min_MgO=params['backtrack']['min_MgO'], 
    min_SiO2=params['backtrack']['min_SiO2'])
if params['backtrack']['apply']:
    s.backtrack_compositions()
if params['PT']['apply']:
    s.compute_pressure_temperature()

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
max_P = s.PT['P'].max() + 1.
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# ---- Fit Tp to suite
if params['suite_Tp']['apply']:
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
        
# ---- Output
if (
    params['backtrack']['save'] or
    params['PT']['save'] or
    params['suite_Tp']['save'] or 
    params['individual_Tp']['save']):
    
    s.write_to_csv(
        params['output_file'],
        write_primary=params['backtrack']['save'],
        write_primary_extended=params['backtrack']['save_extended'],
        write_PT=params['PT']['save'],
        write_suite_Tp=params['suite_Tp']['save'])