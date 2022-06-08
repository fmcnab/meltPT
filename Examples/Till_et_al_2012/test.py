from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- check against plank & forsyth supplementary
s = Suite("TGK12_data.csv", src_FeIII_totFe=0.)

processed = []
for i in range(len(s.data)):
    oxide_wt_hydrous = {}
    for ox in MAJOR_OXIDES:
        oxide_wt_hydrous[ox] = s.data[ox].iloc[i]
    oxide_wt_hydrous = normalise(oxide_wt_hydrous)
    oxide_mole_hydrous = oxide_wt_to_oxide_mole(oxide_wt_hydrous)
    
    oxide_wt_anhydrous = {phase:oxide_wt_hydrous[phase] for phase in oxide_wt_hydrous.keys() if phase != "H2O"}
    oxide_wt_anhydrous = normalise(oxide_wt_anhydrous)
    oxide_mole_anhydrous = {phase:oxide_mole_hydrous[phase] for phase in oxide_mole_hydrous.keys() if phase != "H2O"}
    oxide_mole_anhydrous = normalise(oxide_mole_anhydrous)

    primary_oxide = {}
    for phase in oxide_wt_hydrous:
        primary_oxide[phase + "_primary_wt"] = oxide_wt_hydrous[phase]
    for phase in oxide_wt_anhydrous:
        primary_oxide[phase + "_primary_wt_dry"] = oxide_wt_anhydrous[phase]
    for phase in oxide_mole_hydrous:
        primary_oxide[phase + "_primary_mol"] = oxide_mole_hydrous[phase]
    for phase in oxide_mole_anhydrous:
        primary_oxide[phase + "_primary_mol_dry"] = oxide_mole_anhydrous[phase]

    processed.append(primary_oxide)

s.primary = pd.DataFrame(processed)


def Mg_num(df):
    return df['MgO_primary_mol_dry'] / (df['MgO_primary_mol_dry'] + df['FeO_primary_mol_dry'])

def NaK_num(df):
    return (
        (df['Na2O_primary_mol_dry'] + df['K2O_primary_mol_dry']) / 
        (df['Na2O_primary_mol_dry'] + df['K2O_primary_mol_dry'] + df['CaO_primary_mol_dry'])
        )

def compute_pressure(df, Oliv):
    P = -0.862
    P += 9.471 * Oliv
    P -= 2.383 * (1. - Mg_num(df)) 
    P += 2.922 * NaK_num(df)
    P += 0.218 * df['TiO2_primary_wt_dry'] 
    P -= 0.146 * df['K2O_primary_wt_dry']
    return float(P)

def compute_temperature(df, P):
    T = 1212.
    T += 119.9 * P
    T -= 97.33 * (1. - Mg_num(df)) 
    T -= 87.76 * NaK_num(df)
    T += 3.44 * df['TiO2_primary_wt_dry'] 
    T -= 4.58 * df['K2O_primary_wt_dry']
    return float(T)

def compute_olivine(df, P):
    Oliv = 0.123
    Oliv += 0.085 * P
    Oliv += 0.203 * (1. - Mg_num(df))
    Oliv -= 0.268 * NaK_num(df)
    Oliv -= 0.019 * df['TiO2_primary_wt_dry']
    Oliv += 0.012 * df['K2O_primary_wt_dry']
    return float(Oliv)

def to_minimize(guessP, df):
    # T = compute_temperature(df, P)
    Oliv = compute_olivine(df, guessP)
    newP = compute_pressure(df, Oliv)
    return guessP - newP

from scipy.optimize import root_scalar

fit = root_scalar(to_minimize, args=(s.primary.iloc[0],), x0=1., x1=2.)

P = fit.root
T = compute_temperature(s.primary.iloc[0], P)
