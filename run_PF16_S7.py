from meltPT import *

# ---- check against plank & forsyth supplementary

# # just updating dataframe
# df = parse_csv("PF16_S7.csv", src_FeIII_totFe=0.19)
# df = backtrack_compositions(df, Kd=0.3, verbose=True)
# df = compute_pressure_temperature(df)

# with a class wrapper
s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s.backtrack_compositions(Kd=0.3, verbose=True)
s.compute_pressure_temperature()
